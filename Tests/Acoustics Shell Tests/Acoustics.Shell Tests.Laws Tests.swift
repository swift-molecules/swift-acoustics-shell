import Acoustics_Shell
import Testing

extension `Acoustics.Shell Tests` {
    @Suite struct Laws {}
}

extension `Acoustics.Shell Tests`.Laws {
    typealias T = `Acoustics.Shell Tests`

    /// A quartersawn Sitka spruce panel, 0.40 × 0.25 m × 3 mm, its grain along x.
    static func panel(columns: Int = 16, rows: Int = 10) throws -> (mesh: Mesh<Point<3, Double>>, rim: [Int]) {
        let mesh = try T.rectangle(0.4, 0.25, columns: columns, rows: rows)
        return (mesh, T.rim(mesh, 0.4, 0.25))
    }

    static let spruce = Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter)

    @Test
    func `A flat shell bends as the plate does`() throws {
        let (mesh, rim) = try Self.panel()
        let shell = Acoustics.Shell(
            mesh: mesh,
            thickness: T.metres(0.003),
            material: Self.spruce,
            grain: Vector(x: 1, y: 0, z: 0),
            held: Dictionary(uniqueKeysWithValues: rim.map { ($0, Acoustics.Shell.Hold.clamped) })
        )
        let flat = try Mesh<Point<2, Double>>(
            vertices: mesh.vertices.map { Point(x: $0.x, y: $0.y) },
            faces: mesh.faces.map { face in
                try Polygon(vertices: face.vertices.map { Index(Ordinal(UInt(Int(bitPattern: $0)))) })
            }
        )
        let plate = Acoustics.Plate(mesh: flat, thickness: T.metres(0.003), material: Self.spruce, edge: .clamped)
        let shellModes = try Acoustics.Shell.Modes(shell, count: 6)
        let plateModes = try Acoustics.Plate.Modes(plate, count: 6)
        for (s, p) in zip(shellModes, plateModes) {
            #expect(abs(s.frequency.underlying / p.frequency.underlying - 1) < 1e-6)
            // The shape is the plate's deflection along z, and the modal mass the same at that scale.
            #expect(abs(s.mass.underlying / p.mass.underlying - 1) < 1e-6)
            for (motion, w) in zip(s.shape, p.shape) {
                #expect(abs(abs(motion.z) - abs(w)) < 1e-6)
                #expect(abs(motion.x) + abs(motion.y) < 1e-9)
            }
        }
    }

    @Test
    func `A shell held nowhere has six rigid motions first`() throws {
        let (mesh, _) = try Self.panel(columns: 8, rows: 5)
        let bent = try Mesh(
            vertices: mesh.vertices.map { p in Point(x: p.x, y: p.y, z: 0.3 * (p.x - 0.2) * (p.x - 0.2)) },
            faces: mesh.faces
        )
        let shell = Acoustics.Shell(
            mesh: bent,
            thickness: T.metres(0.003),
            material: Self.spruce,
            grain: Vector(x: 1, y: 0, z: 0)
        )
        let modes = try Acoustics.Shell.Modes(shell, count: 8)
        let f = T.frequencies(modes)
        for value in f.prefix(6) { #expect(value < 1e-3 * f[6], "\(f)") }
        #expect(f[6] > 1)
    }

    @Test
    func `A shell turned in space rings as before`() throws {
        let (mesh, _) = try Self.panel(columns: 8, rows: 5)
        func curved(_ p: Point<3, Double>) -> Vector<3, Double> {
            Vector(x: p.x, y: p.y, z: 0.5 * (p.y - 0.125) * (p.y - 0.125))
        }
        // A turn of 40° about (1, 2, 2) / 3, by Rodrigues' formula.
        let axis = Vector<3, Double>(x: 1.0 / 3, y: 2.0 / 3, z: 2.0 / 3)
        let (c, s) = (0.766044443118978, 0.6427876096865393)
        func turn(_ v: Vector<3, Double>) -> Vector<3, Double> {
            v.scaled(by: c) + axis.cross(v).scaled(by: s) + axis.scaled(by: axis.dot(v) * (1 - c))
        }
        func shell(_ place: (Vector<3, Double>) -> Vector<3, Double>) throws -> Acoustics.Shell {
            let moved = try Mesh(
                vertices: mesh.vertices.map { Point(coordinates: place(curved($0))) },
                faces: mesh.faces
            )
            let corners = [0, 8, mesh.vertices.count - 9, mesh.vertices.count - 1]
            return Acoustics.Shell(
                mesh: moved,
                thickness: T.metres(0.003),
                material: Self.spruce,
                grain: place(Vector(x: 1, y: 0, z: 0)),
                bars: [
                    Acoustics.Shell.Bar(
                        from: Point(coordinates: place(curved(Point(x: 0.17, y: 0, z: 0)))),
                        to: Point(coordinates: place(curved(Point(x: 0.17, y: 0.25, z: 0)))),
                        width: T.metres(0.005),
                        height: T.metres(0.008),
                        material: Self.spruce,
                        side: place(Vector(x: 0, y: 0, z: -1))
                    )
                ],
                held: Dictionary(uniqueKeysWithValues: corners.map { ($0, Acoustics.Shell.Hold.pinned) })
            )
        }
        let still = T.frequencies(try Acoustics.Shell.Modes(try shell { $0 }, count: 6))
        let turned = T.frequencies(try Acoustics.Shell.Modes(try shell(turn), count: 6))
        for (a, b) in zip(still, turned) { #expect(abs(a / b - 1) < 1e-8, "\(a) \(b)") }
    }

    @Test
    func `A shell scaled in every length rings at the scale's inverse`() throws {
        let (mesh, rim) = try Self.panel(columns: 8, rows: 5)
        func shell(_ scale: Double) throws -> Acoustics.Shell {
            let curved = try Mesh(
                vertices: mesh.vertices.map { p in
                    Point(x: scale * p.x, y: scale * p.y, z: scale * 0.4 * (p.x - 0.2) * (p.x - 0.2))
                },
                faces: mesh.faces
            )
            return Acoustics.Shell(
                mesh: curved,
                thickness: T.metres(0.003 * scale),
                material: Self.spruce,
                grain: Vector(x: 1, y: 0, z: 0),
                held: Dictionary(uniqueKeysWithValues: rim.map { ($0, Acoustics.Shell.Hold.pinned) })
            )
        }
        let one = try Acoustics.Shell.Modes(try shell(1), count: 4)
        let two = try Acoustics.Shell.Modes(try shell(2), count: 4)
        for (a, b) in zip(one, two) {
            #expect(abs(b.frequency.underlying * 2 / a.frequency.underlying - 1) < 1e-8)
            #expect(abs(b.mass.underlying / (8 * a.mass.underlying) - 1) < 1e-8)
        }
    }

    @Test
    func `Bad input is refused`() throws {
        let (mesh, _) = try Self.panel(columns: 2, rows: 2)
        let shell = Acoustics.Shell(
            mesh: mesh,
            thickness: T.metres(0.003),
            material: Self.spruce,
            grain: Vector(x: 1, y: 0, z: 0)
        )
        var stray = shell
        stray.materials = [9: Self.spruce]
        #expect(throws: Acoustics.Shell.Error.face(9)) { try Acoustics.Shell.Modes(stray, count: 1) }
        var held = shell
        held.held = [99: .clamped]
        #expect(throws: Acoustics.Shell.Error.vertex(99)) { try Acoustics.Shell.Modes(held, count: 1) }
        var bar = shell
        bar.bars = [
            Acoustics.Shell.Bar(
                from: Point(x: 0.1, y: 0, z: 0),
                to: Point(x: 0.1, y: 0, z: 0),
                width: T.metres(0.005),
                height: T.metres(0.008),
                material: Self.spruce,
                side: Vector(x: 0, y: 0, z: 1)
            )
        ]
        #expect(throws: Acoustics.Shell.Error.bar(0)) { try Acoustics.Shell.Modes(bar, count: 1) }
        #expect(throws: Acoustics.Shell.Error.count) { try Acoustics.Shell.Modes(shell, count: 1_000) }
    }
}
