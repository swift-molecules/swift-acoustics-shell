public import Acoustics_Plate
internal import Galerkin
internal import Sparse

extension Acoustics.Shell {
    /// A shell's lowest modes, computed once from its shape and material and kept as a value.
    public struct Modes {
        /// The shell they belong to.
        public let shell: Acoustics.Shell
        /// The modes, lowest first; a repeated frequency appears once for each mode that has it.
        public let elements: [Mode]

        /// The `count` lowest modes of `shell`.
        ///
        /// A shell held nowhere has six modes at zero frequency to rounding first: its rigid motions.
        public init(_ shell: Acoustics.Shell, count: Int) throws(Acoustics.Shell.Error) {
            let system = try Acoustics.Shell.System(shell)
            guard count >= 0, count <= system.order else { throw .count }
            let pairs: [Eigen.Pair<Double>]
            do throws(Eigen.Error) {
                let pencil = try Eigen.Pencil(system.stiffness, system.mass)
                pairs = try pencil.lowest(Cardinal(UInt(count)))
            } catch {
                throw .eigen(error)
            }
            self.shell = shell
            self.elements = pairs.map { system.mode($0) }
        }
    }
}

extension Acoustics.Shell.Modes: Hashable {}
extension Acoustics.Shell.Modes: Sendable {}

extension Acoustics.Shell.Modes: RandomAccessCollection {
    public var startIndex: Int { elements.startIndex }
    public var endIndex: Int { elements.endIndex }
    public subscript(position: Int) -> Acoustics.Shell.Mode { elements[position] }
}

extension Acoustics.Shell {
    /// A shell's stiffness and mass over its unknowns that are not held, and what turns an eigenvector back into a
    /// mode.
    internal struct System {
        let stiffness: Sparse.Symmetric<Double>
        let mass: Sparse.Symmetric<Double>
        /// For each vertex, the unknowns of its motion along x, y and z, or nil where held.
        let motion: [[Int?]]

        var order: Int { stiffness.order }
    }
}

extension Acoustics.Shell.System {
    init(_ shell: Acoustics.Shell) throws(Acoustics.Shell.Error) {
        let points = shell.mesh.vertices
        var facets: [Acoustics.Shell.Facet] = []
        for (index, face) in shell.mesh.faces.enumerated() {
            let corners = face.vertices.map { Int(bitPattern: $0) }
            guard corners.count == 4,
                let facet = Acoustics.Shell.Facet(corners: corners, points: points, grain: shell.grain)
            else { throw .face(index) }
            facets.append(facet)
        }
        if let stray = shell.materials.keys.first(where: { !facets.indices.contains($0) }) { throw .face(stray) }
        if let stray = shell.thicknesses.keys.first(where: { !facets.indices.contains($0) }) { throw .face(stray) }
        if let stray = shell.held.keys.first(where: { !points.indices.contains($0) }) { throw .vertex(stray) }
        for (index, bar) in shell.bars.enumerated() {
            let length = (bar.end.coordinates - bar.start.coordinates).length
            guard length > 0, bar.width.underlying > 0, bar.height.underlying > 0, bar.side.length > 0 else {
                throw .bar(index)
            }
        }
        var used = [Bool](repeating: false, count: points.count)
        for facet in facets {
            for corner in facet.corners { used[corner] = true }
        }
        var held: [Int: Set<Int>] = [:]
        for (vertex, hold) in shell.held { held[vertex] = hold.components }
        for vertex in points.indices where !used[vertex] { held[vertex] = Set(0..<6) }
        let unknowns = Galerkin.Unknowns(vertices: points.count, components: 6, held: held)
        var stiffness = Galerkin.Assembly<Double>(unknowns)
        var mass = stiffness
        for (index, facet) in facets.enumerated() {
            let h = (shell.thicknesses[index] ?? shell.thickness).underlying
            var (k, m) = facet.matrices(thickness: h, material: shell.materials[index] ?? shell.material)
            for bar in shell.bars {
                guard let (kb, mb) = facet.bar(bar, thickness: h) else { continue }
                for entry in 0..<576 {
                    k[entry] += kb[entry]
                    m[entry] += mb[entry]
                }
            }
            stiffness.add(facet.turned(k), at: facet.corners)
            mass.add(facet.turned(m), at: facet.corners)
        }
        self.stiffness = stiffness.matrix
        self.mass = mass.matrix
        self.motion = points.indices.map { vertex in (0..<3).map { unknowns[vertex, $0] } }
    }

    /// The mode of an eigenpair: its frequency, its shape scaled so that the longest motion is one, and its mass at
    /// that scale.
    func mode(_ pair: Eigen.Pair<Double>) -> Acoustics.Shell.Mode {
        var shape = motion.map { components in
            let values = components.map { $0.map { pair.vector[$0] } ?? 0 }
            return Vector(x: values[0], y: values[1], z: values[2])
        }
        var peak = 0.0
        var sign = 1.0
        var largest = 0.0
        for motion in shape {
            peak = Swift.max(peak, motion.length)
            for component in [motion.x, motion.y, motion.z] where Swift.abs(component) > largest {
                largest = Swift.abs(component)
                sign = component < 0 ? -1 : 1
            }
        }
        if peak > 0 {
            shape = shape.map { $0.scaled(by: sign / peak) }
        }
        let frequency = Swift.max(pair.value, 0).squareRoot() / (2 * Double.pi)
        // The vector is B-normalized, so at a longest motion of one the modal mass is 1 / peak².
        let mass = peak > 0 ? 1 / (peak * peak) : 0
        return Acoustics.Shell.Mode(frequency: Frequency.hertz(frequency), shape: shape, mass: Mass.kilograms(mass))
    }
}

extension Acoustics.Shell {
    /// The unit normal at each vertex: the mean of the normals of the faces that use it, each weighted by the face's
    /// area and pointing to where its corners turn counterclockwise. A vertex no face uses has zero.
    public var normals: [Vector<3, Double>] {
        let points = mesh.vertices
        var sums = [Vector<3, Double>](repeating: .zero, count: points.count)
        let middle = Galerkin.Reference.Quadrilateral.shape(at: Point(x: 0.0, y: 0.0))
        for face in mesh.faces {
            let corners = face.vertices.map { Int(bitPattern: $0) }
            guard corners.count == 4 else { continue }
            let normal = Galerkin.Jacobian.Surface(middle, corners: corners.map { points[$0] }).normal
            for corner in corners { sums[corner] += normal }
        }
        return sums.map { $0.length > 0 ? $0.scaled(by: 1 / $0.length) : $0 }
    }
}
