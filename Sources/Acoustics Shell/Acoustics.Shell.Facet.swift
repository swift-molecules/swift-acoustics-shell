internal import Galerkin

extension Acoustics.Shell {
    /// One flat facet of a shell: a quadrilateral laid into the plane through its centre square to its normal there,
    /// with axes e₁ along the grain, e₂ across it and e₃ along the normal.
    ///
    /// Each corner has six unknowns in space, the motion (u, v, w) and the turn (θx, θy, θz). In the facet's axes
    /// the bending is the plate's MITC4 element, with βx = θ₂ and βy = −θ₁; the membrane is the bilinear plane-stress
    /// quadrilateral at 2 × 2 Gauss points; and the turn about the normal, which neither resists, is tied to the
    /// membrane's own rotation ½(∂v/∂x − ∂u/∂y) by a weak penalty (Hughes and Brezzi 1989), so that a facet that lies
    /// flat among its neighbours keeps no free turn.
    internal struct Facet {
        /// The corners, as the mesh's vertices, in the order whose normal is e₃.
        let corners: [Int]
        /// The facet's axes in space.
        let axes: (first: Vector<3, Double>, second: Vector<3, Double>, normal: Vector<3, Double>)
        /// The centre the corners are measured from.
        let centre: Vector<3, Double>
        /// The plate element over the corners in the facet's axes.
        let plate: Acoustics.Plate.Element
    }
}

extension Acoustics.Shell.Facet {
    typealias Square = Galerkin.Reference.Quadrilateral

    static let quadrature = Square.quadrature(degree: 3, Double.self)

    /// The share of the in-plane shear stiffness that ties the turn about the normal to the membrane's rotation.
    static let drilling = 1e-3

    /// The facet over these corners, or nil where it has no area or turns over when laid flat.
    init?(corners: [Int], points: [Point<3, Double>], grain: Vector<3, Double>) {
        let positions = corners.map { points[$0] }
        let middle = Square.shape(at: Point(x: 0.0, y: 0.0))
        let jacobian = Galerkin.Jacobian.Surface(middle, corners: positions)
        let length = jacobian.normal.length
        guard length > 0 else { return nil }
        let normal = jacobian.normal.scaled(by: 1 / length)
        var first = grain - normal.scaled(by: grain.dot(normal))
        if first.length < 1e-6 * grain.length || grain.length == 0 {
            first = jacobian.xi - normal.scaled(by: jacobian.xi.dot(normal))
        }
        first = first.scaled(by: 1 / first.length)
        let second = normal.cross(first)
        var centre = Vector<3, Double>.zero
        for position in positions { centre += position.coordinates.scaled(by: 0.25) }
        let flat = positions.map { position in
            let offset = position.coordinates - centre
            return Point(x: offset.dot(first), y: offset.dot(second))
        }
        let plate = Acoustics.Plate.Element(flat)
        guard plate.determinant > 0 else { return nil }
        self.corners = corners
        self.axes = (first, second, normal)
        self.centre = centre
        self.plate = plate
    }

    /// A point of space in the facet's axes: along e₁, along e₂, and off the facet along e₃.
    func local(_ point: Point<3, Double>) -> (x: Double, y: Double, z: Double) {
        let offset = point.coordinates - centre
        return (offset.dot(axes.first), offset.dot(axes.second), offset.dot(axes.normal))
    }

    /// The shape functions and their gradients in the facet's axes at (ξ, η), with the area element.
    func gradients(at point: Point<2, Double>) -> (values: [Double], x: [Double], y: [Double], area: Double) {
        let shape = Square.shape(at: point)
        let jacobian = Galerkin.Jacobian(shape, corners: plate.corners)
        let (x, y) = jacobian.gradients(shape)
        return (shape.values, x, y, jacobian.determinant)
    }

    /// The reference point (ξ, η) that maps onto (x, y) in the facet's plane, by Newton's method from the centre.
    func reference(x: Double, y: Double) -> Point<2, Double> {
        var xi = 0.0
        var eta = 0.0
        for _ in 0..<20 {
            let shape = Square.shape(at: Point(x: xi, y: eta))
            let jacobian = Galerkin.Jacobian(shape, corners: plate.corners)
            var px = 0.0
            var py = 0.0
            for (value, corner) in zip(shape.values, plate.corners) {
                px += value * corner.x
                py += value * corner.y
            }
            let (dx, dy) = (x - px, y - py)
            let det = jacobian.determinant
            let dxi = (jacobian.yEta * dx - jacobian.xEta * dy) / det
            let deta = (-jacobian.yXi * dx + jacobian.xXi * dy) / det
            xi += dxi
            eta += deta
            if Swift.abs(dxi) + Swift.abs(deta) < 1e-13 { break }
        }
        return Point(x: xi, y: eta)
    }
}

extension Acoustics.Shell.Facet {
    /// The 24 × 24 stiffness and mass, row-major over (u, v, w, θ₁, θ₂, θ₃) at each corner in the facet's axes.
    func matrices(
        thickness h: Double,
        material: Acoustics.Plate.Material
    ) -> (stiffness: [Double], mass: [Double]) {
        var k = [Double](repeating: 0, count: 576)
        var m = [Double](repeating: 0, count: 576)
        let (bending, inertia) = plate.matrices(thickness: Space.Physical.metres(h), material: material)
        // The plate's (w, βx, βy) are (w, θ₂, −θ₁).
        let place = [(2, 1.0), (4, 1.0), (3, -1.0)]
        for a in 0..<12 {
            let (ra, sa) = place[a % 3]
            let row = 6 * (a / 3) + ra
            for b in 0..<12 {
                let (rb, sb) = place[b % 3]
                let column = 6 * (b / 3) + rb
                k[24 * row + column] += sa * sb * bending[12 * a + b]
                m[24 * row + column] += sa * sb * inertia[12 * a + b]
            }
        }
        let ex = material.young.x.underlying
        let ey = material.young.y.underlying
        let nxy = material.poisson.underlying
        let denominator = 1 - nxy * nxy * ey / ex
        let q11 = h * ex / denominator
        let q12 = h * nxy * ey / denominator
        let q22 = h * ey / denominator
        let q66 = h * material.shear.xy.underlying
        let penalty = Self.drilling * q66
        let rho = material.density.underlying
        for (point, weight) in zip(Self.quadrature.points, Self.quadrature.weights) {
            let g = gradients(at: point)
            let dA = weight * g.area
            var ex = [Double](repeating: 0, count: 24)
            var ey = ex
            var gxy = ex
            var drill = ex
            for i in 0..<4 {
                ex[6 * i] = g.x[i]
                ey[6 * i + 1] = g.y[i]
                gxy[6 * i] = g.y[i]
                gxy[6 * i + 1] = g.x[i]
                drill[6 * i + 5] = g.values[i]
                drill[6 * i] = 0.5 * g.y[i]
                drill[6 * i + 1] = -0.5 * g.x[i]
            }
            for a in 0..<24 where ex[a] != 0 || ey[a] != 0 || gxy[a] != 0 || drill[a] != 0 {
                let nx = q11 * ex[a] + q12 * ey[a]
                let ny = q12 * ex[a] + q22 * ey[a]
                let nxy = q66 * gxy[a]
                let p = penalty * drill[a]
                for b in 0..<24 {
                    k[24 * a + b] += dA * (nx * ex[b] + ny * ey[b] + nxy * gxy[b] + p * drill[b])
                }
            }
            for i in 0..<4 {
                for j in 0..<4 {
                    let nn = dA * g.values[i] * g.values[j]
                    m[24 * (6 * i) + 6 * j] += rho * h * nn
                    m[24 * (6 * i + 1) + 6 * j + 1] += rho * h * nn
                    m[24 * (6 * i + 5) + 6 * j + 5] += rho * h * h * h / 12 * nn
                }
            }
        }
        return (k, m)
    }
}

extension Acoustics.Shell.Facet {
    /// The part of the segment from `start` to `end`, laid into the facet's plane, that lies inside it: its ends
    /// in the facet's axes and the share of it that counts here, a half where it runs along the facet's edge.
    func clip(_ start: Point<3, Double>, _ end: Point<3, Double>) -> (from: (Double, Double), to: (Double, Double), share: Double)? {
        let a = local(start)
        let b = local(end)
        var low = 0.0
        var high = 1.0
        let corners = plate.corners
        let size = (corners[2].coordinates - corners[0].coordinates).length
        let tolerance = 1e-9 * size
        var along = false
        for index in 0..<4 {
            let p = corners[index]
            let q = corners[(index + 1) % 4]
            // The inward normal of the edge from p to q of a counterclockwise polygon is (−(q.y − p.y), q.x − p.x).
            let nx = -(q.y - p.y)
            let ny = q.x - p.x
            let length = (nx * nx + ny * ny).squareRoot()
            let start = (nx * (a.x - p.x) + ny * (a.y - p.y)) / length
            let rate = (nx * (b.x - a.x) + ny * (b.y - a.y)) / length
            if Swift.abs(rate) < tolerance {
                if start < -tolerance { return nil }
                if start < tolerance { along = true }
                continue
            }
            let t = -start / rate
            if rate > 0 { low = Swift.max(low, t) } else { high = Swift.min(high, t) }
        }
        guard high - low > 1e-12 else { return nil }
        let from = (a.x + low * (b.x - a.x), a.y + low * (b.y - a.y))
        let to = (a.x + high * (b.x - a.x), a.y + high * (b.y - a.y))
        return (from, to, along ? 0.5 : 1)
    }

    /// The 24 × 24 stiffness and mass a bar adds over the part of it inside the facet, in the facet's axes, for a
    /// shell of thickness `h` here.
    func bar(_ bar: Acoustics.Shell.Bar, thickness h: Double) -> (stiffness: [Double], mass: [Double])? {
        guard let piece = clip(bar.start, bar.end) else { return nil }
        let dx = piece.to.0 - piece.from.0
        let dy = piece.to.1 - piece.from.1
        let length = (dx * dx + dy * dy).squareRoot()
        guard length > 0 else { return nil }
        let (sx, sy) = (dx / length, dy / length)
        let b = bar.width.underlying
        let t = bar.height.underlying
        let sign: Double = bar.side.dot(axes.normal) < 0 ? -1 : 1
        let e = sign * (h + t) / 2
        let area = b * t
        let own = b * t * t * t / 12
        let (thin, deep) = (Swift.min(b, t), Swift.max(b, t))
        let r = thin / deep
        let torsion = deep * thin * thin * thin * (1.0 / 3 - 0.21 * r * (1 - r * r * r * r / 12))
        let young = bar.material.young.x.underlying
        let shear = bar.material.shear.xz.underlying
        let twist = bar.material.shear.xy.underlying
        let rho = bar.material.density.underlying
        var k = [Double](repeating: 0, count: 576)
        var m = [Double](repeating: 0, count: 576)
        func add(_ matrix: inout [Double], _ row: [Double], _ scale: Double) {
            for a in 0..<24 where row[a] != 0 {
                for c in 0..<24 where row[c] != 0 { matrix[24 * a + c] += scale * row[a] * row[c] }
            }
        }
        func rows(at s: Double) -> (g: (values: [Double], x: [Double], y: [Double], area: Double), axial: [Double], curvature: [Double], twist: [Double], shear: [Double]) {
            let point = reference(x: piece.from.0 + s * dx, y: piece.from.1 + s * dy)
            let g = gradients(at: point)
            var axial = [Double](repeating: 0, count: 24)
            var curvature = axial
            var turning = axial
            var shearing = axial
            for i in 0..<4 {
                let along = sx * g.x[i] + sy * g.y[i]
                // ε_s = s·∇(u·s); κ_s = s·∇(β·s) with βx = θ₂, βy = −θ₁; τ = s·∇(θ·s); γ = w,s + β·s.
                axial[6 * i] = sx * along
                axial[6 * i + 1] = sy * along
                curvature[6 * i + 4] = sx * along
                curvature[6 * i + 3] = -sy * along
                turning[6 * i + 3] = sx * along
                turning[6 * i + 4] = sy * along
                shearing[6 * i + 2] = along
                shearing[6 * i + 4] = g.values[i] * sx
                shearing[6 * i + 3] = -g.values[i] * sy
            }
            return (g, axial, curvature, turning, shearing)
        }
        let rule = Galerkin.Quadrature<1, Double>.gauss(count: 3)
        for (node, weight) in zip(rule.nodes, rule.weights) {
            let s = (node + 1) / 2
            let ds = piece.share * weight * length / 2
            let at = rows(at: s)
            var stretch = at.axial
            for a in 0..<24 { stretch[a] += e * at.curvature[a] }
            add(&k, stretch, ds * young * area)
            add(&k, at.curvature, ds * young * own)
            add(&k, at.twist, ds * twist * torsion)
            var mx = [Double](repeating: 0, count: 24)
            var my = mx
            var mz = mx
            for i in 0..<4 {
                let n = at.g.values[i]
                mx[6 * i] = n
                mx[6 * i + 4] = e * n
                my[6 * i + 1] = n
                my[6 * i + 3] = -e * n
                mz[6 * i + 2] = n
            }
            add(&m, mx, ds * rho * area)
            add(&m, my, ds * rho * area)
            add(&m, mz, ds * rho * area)
        }
        // The bar's shear at the middle of the piece alone, so that a short piece does not lock.
        let middle = rows(at: 0.5)
        add(&k, middle.shear, piece.share * length * 5.0 / 6 * shear * area)
        return (k, m)
    }
}

extension Acoustics.Shell.Facet {
    /// `matrix` over the facet's axes turned into one over the axes of space: Tᵀ A T, T taking each corner's motion
    /// and turn in space to the facet's axes.
    func turned(_ matrix: [Double]) -> [Double] {
        let rows = [axes.first, axes.second, axes.normal]
        // R[i][j] = eᵢ · xⱼ
        var r = [Double](repeating: 0, count: 9)
        for i in 0..<3 {
            r[3 * i] = rows[i].x
            r[3 * i + 1] = rows[i].y
            r[3 * i + 2] = rows[i].z
        }
        // A T, block by block of three.
        var at = [Double](repeating: 0, count: 576)
        for row in 0..<24 {
            for block in 0..<8 {
                for j in 0..<3 {
                    var sum = 0.0
                    for i in 0..<3 { sum += matrix[24 * row + 3 * block + i] * r[3 * i + j] }
                    at[24 * row + 3 * block + j] = sum
                }
            }
        }
        var result = [Double](repeating: 0, count: 576)
        for block in 0..<8 {
            for j in 0..<3 {
                for column in 0..<24 {
                    var sum = 0.0
                    for i in 0..<3 { sum += r[3 * i + j] * at[24 * (3 * block + i) + column] }
                    result[24 * (3 * block + j) + column] = sum
                }
            }
        }
        return result
    }
}
