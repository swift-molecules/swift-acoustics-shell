import Acoustics_Shell
import Numeric
import Testing

extension `Acoustics.Shell Tests` {
    /// Shells and barred strips whose modes are known in closed form.
    @Suite struct Closed {}
}

extension `Acoustics.Shell Tests`.Closed {
    typealias T = `Acoustics.Shell Tests`

    /// The eigenvalues of a symmetric 3 × 3 matrix, ascending, by Jacobi's rotations.
    static func eigenvalues(_ matrix: [[Double]]) -> [Double] {
        var a = matrix
        for _ in 0..<100 {
            var off = 0.0
            for p in 0..<3 {
                for q in (p + 1)..<3 { off += a[p][q] * a[p][q] }
            }
            if off < 1e-30 * (a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2]) { break }
            for p in 0..<3 {
                for q in (p + 1)..<3 where a[p][q] != 0 {
                    let theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                    let t = (theta >= 0 ? 1.0 : -1.0) / (abs(theta) + (theta * theta + 1).squareRoot())
                    let c = 1 / (t * t + 1).squareRoot()
                    let s = t * c
                    for k in 0..<3 {
                        let (akp, akq) = (a[k][p], a[k][q])
                        a[k][p] = c * akp - s * akq
                        a[k][q] = s * akp + c * akq
                    }
                    for k in 0..<3 {
                        let (apk, aqk) = (a[p][k], a[q][k])
                        a[p][k] = c * apk - s * aqk
                        a[q][k] = s * apk + c * aqk
                    }
                }
            }
        }
        return [a[0][0], a[1][1], a[2][2]].sorted()
    }

    /// The lowest frequency of a cylinder of radius `a`, length `l` and thickness `h`, its ends held by shear
    /// diaphragms (v = w = 0), in the mode with `m` half waves along it and `n` waves round it, by Love's equations.
    ///
    /// Love's strains and changes of curvature (Leissa 1973, Vibration of Shells, NASA SP-288, §2.2; Soedel 2004,
    /// Vibrations of Shells and Plates, 3rd ed., §8.4): ε_x = u,x; ε_θ = (v,θ + w)/a; γ = v,x + u,θ/a;
    /// κ_x = −w,xx; κ_θ = (v,θ − w,θθ)/a²; τ = (v,x − w,xθ)/a. With u = A cos λx cos nθ, v = B sin λx sin nθ and
    /// w = C sin λx cos nθ, λ = mπ/l, the strain and kinetic energies are quadratic forms in (A, B, C) and ω² the
    /// eigenvalues of the first over ρh.
    static func cylinder(
        radius a: Double,
        length l: Double,
        thickness h: Double,
        material: Acoustics.Plate.Material,
        m: Int,
        n: Int
    ) -> Double {
        let e = material.young.x.underlying
        let nu = material.poisson.underlying
        let rho = material.density.underlying
        let membrane = e * h / (1 - nu * nu)
        let bending = e * h * h * h / (12 * (1 - nu * nu))
        let lambda = Double(m) * Double.pi / l
        let nn = Double(n)
        let ex = [-lambda, 0, 0]
        let et = [0, nn / a, 1 / a]
        let g = [-nn / a, lambda, 0]
        let kx = [0, 0, lambda * lambda]
        let kt = [0, nn / (a * a), nn * nn / (a * a)]
        let tw = [0, lambda / a, lambda * nn / a]
        var k = [[Double]](repeating: [0, 0, 0], count: 3)
        for i in 0..<3 {
            for j in 0..<3 {
                let normal = ex[i] * ex[j] + et[i] * et[j] + nu * (ex[i] * et[j] + et[i] * ex[j])
                let stretch = normal + (1 - nu) / 2 * g[i] * g[j]
                let curving = kx[i] * kx[j] + kt[i] * kt[j] + nu * (kx[i] * kt[j] + kt[i] * kx[j])
                let bend = curving + 2 * (1 - nu) * tw[i] * tw[j]
                k[i][j] = membrane * stretch + bending * bend
            }
        }
        let omega2 = eigenvalues(k)[0] / (rho * h)
        return omega2.squareRoot() / (2 * Double.pi)
    }

    /// The modes of the cylinder of radius 0.1 m, length 0.2 m and thickness 1 mm, `around` × `along` facets, each
    /// with the half waves m along it and the waves n round it that its shape shows, and its frequency.
    static func cylinder(around: Int, along: Int) throws -> [(m: Int, n: Int, frequency: Double)] {
        let (a, l, h) = (0.1, 0.2, 0.001)
        let mesh = try T.grid(columns: around, rows: along, wrap: true) { i, j in
            let theta = 2 * Double.pi * Double(i) / Double(around)
            return Point(x: a * Double.cos(theta), y: a * Double.sin(theta), z: l * Double(j) / Double(along))
        }
        // The ends hold the motion across the axis and leave the axial motion and every turn free.
        let ends = mesh.vertices.indices.filter {
            abs(mesh.vertices[$0].z) < 1e-12 || abs(mesh.vertices[$0].z - l) < 1e-12
        }
        let shell = Acoustics.Shell(
            mesh: mesh,
            thickness: T.metres(h),
            material: T.steel,
            grain: Vector(x: 0, y: 0, z: 1),
            held: Dictionary(uniqueKeysWithValues: ends.map { ($0, [.x, .y] as Acoustics.Shell.Hold) })
        )
        let modes = try Acoustics.Shell.Modes(shell, count: 17)
        // The first is the free sliding along the axis.
        #expect(modes[0].frequency.underlying < 1e-3 * modes[1].frequency.underlying)
        func radial(_ mode: Acoustics.Shell.Mode, _ i: Int, _ j: Int) -> Double {
            let theta = 2 * Double.pi * Double(i % around) / Double(around)
            return mode.shape[j * around + i % around].dot(Vector(x: Double.cos(theta), y: Double.sin(theta), z: 0))
        }
        func changes(_ values: [Double]) -> Int {
            let peak = values.map { abs($0) }.max() ?? 0
            let signs = values.filter { abs($0) > 1e-3 * peak }.map { $0 > 0 }
            return zip(signs, signs.dropFirst()).filter { $0 != $1 }.count
        }
        return modes.dropFirst().map { mode in
            // Round the ring a quarter of the way along, where no mode with m ≤ 3 has a node.
            let ring = (0...around).map { radial(mode, $0, along / 4) }
            let n = changes(ring) / 2
            let widest = (0..<around).max { abs(radial(mode, $0, along / 4)) < abs(radial(mode, $1, along / 4)) }!
            let m = changes((1..<along).map { radial(mode, widest, $0) }) + 1
            return (m, n, mode.frequency.underlying)
        }
    }

    @Test
    func `A cylinder held by shear diaphragms rings as Love's equations give`() throws {
        // Every mode of the fine mesh with at least 16 facets to a wave round (n ≤ 6) within 2% of Love's, and every
        // mode's error at least three times smaller than on a mesh twice as coarse each way: second-order convergence
        // onto Love's frequencies. At n = 7 and 8 the fine mesh is still 2–3% high.
        let coarse = try Self.cylinder(around: 48, along: 16)
        let fine = try Self.cylinder(around: 96, along: 32)
        for mode in fine {
            let love = Self.cylinder(
                radius: 0.1,
                length: 0.2,
                thickness: 0.001,
                material: T.steel,
                m: mode.m,
                n: mode.n
            )
            let error = mode.frequency / love - 1
            #expect(mode.n > 6 || abs(error) < 0.02, "m \(mode.m) n \(mode.n): \(mode.frequency) against \(love)")
            if let rough = coarse.first(where: { $0.m == mode.m && $0.n == mode.n }) {
                let before = rough.frequency / love - 1
                #expect(abs(before) > 3 * abs(error), "m \(mode.m) n \(mode.n): \(before) then \(error)")
            }
        }
    }

    /// A sphere of radius `radius` made of six grids of `cells` × `cells` quadrilaterals, the faces of a cube carried
    /// out onto it, every face's corners turning counterclockwise seen from outside.
    static func sphere(radius: Double, cells: Int) throws -> Mesh<Point<3, Double>> {
        var vertices: [Point<3, Double>] = []
        var numbers: [[Int]: Int] = [:]
        var faces: [Polygon<Index<Point<3, Double>>>] = []
        func vertex(_ cube: [Int]) -> Int {
            if let number = numbers[cube] { return number }
            let v = Vector<3, Double>(
                x: Double(cube[0]) / Double(cells) - 1,
                y: Double(cube[1]) / Double(cells) - 1,
                z: Double(cube[2]) / Double(cells) - 1
            )
            vertices.append(Point(coordinates: v.scaled(by: radius / v.length)))
            numbers[cube] = vertices.count - 1
            return vertices.count - 1
        }
        let top = 2 * cells
        for axis in 0..<3 {
            for level in [0, top] {
                for i in 0..<top {
                    for j in 0..<top {
                        func at(_ p: Int, _ q: Int) -> Int {
                            var cube = [0, 0, 0]
                            cube[axis] = level
                            cube[(axis + 1) % 3] = p
                            cube[(axis + 2) % 3] = q
                            return vertex(cube)
                        }
                        var corners = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]
                        if level == 0 { corners.reverse() }
                        faces.append(try Polygon(vertices: corners.map { Index(Ordinal(UInt($0))) }))
                    }
                }
            }
        }
        return try Mesh(vertices: vertices, faces: faces)
    }

    @Test
    func `A thin sphere rings as Lamb's membrane equation gives`() throws {
        // Lamb (1882), Proc. London Math. Soc. 14, 50; in the membrane limit, as Kalnins (1964), J. Acoust. Soc. Am.
        // 36, 74, and Soedel (2004) give it: Ω⁴ − (1 + 3ν + λ)Ω² + (λ − 2)(1 − ν²) = 0, λ = n(n + 1),
        // Ω² = ρ(1 − ν²)a²ω²/E, each root 2n + 1 times over. The lower root of n = 2 and n = 3, at a/h = 100, where
        // bending adds less than 0.2%.
        let (a, h) = (0.1, 0.001)
        let mesh = try Self.sphere(radius: a, cells: 6)
        let shell = Acoustics.Shell(
            mesh: mesh,
            thickness: T.metres(h),
            material: T.steel,
            grain: Vector(x: 0, y: 0, z: 1)
        )
        let modes = try Acoustics.Shell.Modes(shell, count: 18)
        let found = Array(T.frequencies(modes).dropFirst(6))
        let e = T.steel.young.x.underlying
        let nu = T.steel.poisson.underlying
        let rho = T.steel.density.underlying
        func lamb(_ n: Int) -> Double {
            let l = Double(n * (n + 1))
            let b = 1 + 3 * nu + l
            let c = (l - 2) * (1 - nu * nu)
            let omega2 = (b - (b * b - 4 * c).squareRoot()) / 2
            return (omega2 * e / (rho * (1 - nu * nu) * a * a)).squareRoot() / (2 * Double.pi)
        }
        // The sphere the facets enclose is smaller than the one they were drawn on; compare with the sphere of equal
        // area, as the membrane frequencies scale with 1 / a.
        var area = 0.0
        for face in mesh.faces {
            let p = face.vertices.map { mesh.vertices[Int(bitPattern: $0)].coordinates }
            area += 0.5 * (p[2] - p[0]).cross(p[3] - p[1]).length
        }
        let equal = (area / (4 * Double.pi)).squareRoot()
        let expected = Array(repeating: lamb(2) * a / equal, count: 5) + Array(repeating: lamb(3) * a / equal, count: 7)
        for (f, x) in zip(found, expected) { #expect(abs(f / x - 1) < 0.02, "\(f) against \(x)") }
    }

    @Test
    func `A barred strip bends as a beam of its whole section`() throws {
        // A strip 0.5 m × 20 mm × 2 mm with a bar 6 mm wide and 12 mm high along its middle line, which is an edge
        // between two faces, simply supported at its ends and held from twisting and moving sideways: a beam whose
        // section is the strip and the bar together, bending about their common centroid, so that
        // fₙ = (n²π / 2L²) √(EI / m′).
        let (l, w, h, b, t) = (0.5, 0.02, 0.002, 0.006, 0.012)
        let wood = T.isotropic(density: 500, young: 10e9, poisson: 0)
        let mesh = try T.rectangle(l, w, columns: 40, rows: 2)
        var held: [Int: Acoustics.Shell.Hold] = [:]
        for (index, p) in mesh.vertices.enumerated() {
            var hold: Acoustics.Shell.Hold = [.y, .turnX]
            if abs(p.x) < 1e-12 { hold.formUnion([.x, .z]) }
            if abs(p.x - l) < 1e-12 { hold.formUnion(.z) }
            held[index] = hold
        }
        let shell = Acoustics.Shell(
            mesh: mesh,
            thickness: T.metres(h),
            material: wood,
            grain: Vector(x: 1, y: 0, z: 0),
            bars: [
                Acoustics.Shell.Bar(
                    from: Point(x: 0, y: w / 2, z: 0),
                    to: Point(x: l, y: w / 2, z: 0),
                    width: T.metres(b),
                    height: T.metres(t),
                    material: wood,
                    side: Vector(x: 0, y: 0, z: 1)
                )
            ],
            held: held
        )
        let modes = try Acoustics.Shell.Modes(shell, count: 2)
        let (a1, a2) = (w * h, b * t)
        let e = (h + t) / 2
        let centroid = a2 * e / (a1 + a2)
        let inertia =
            w * h * h * h / 12 + a1 * centroid * centroid + b * t * t * t / 12 + a2 * (e - centroid) * (e - centroid)
        let perLength = 500 * (a1 + a2)
        for (index, mode) in modes.enumerated() {
            let n = Double(index + 1)
            let expected = n * n * Double.pi / (2 * l * l) * (10e9 * inertia / perLength).squareRoot()
            let found = mode.frequency.underlying
            #expect(abs(found / expected - 1) < 0.01, "\(found) against \(expected)")
        }
    }
}
