import Testing

import Acoustics_Shell

@Suite struct `Acoustics.Shell Tests` {}

extension `Acoustics.Shell Tests` {
    typealias Space3 = Point<3, Double>

    /// Quadrilaterals over a grid of `columns` × `rows` cells, vertex (i, j) at `point(i, j)`. With `wrap`, the last
    /// column of vertices is the first, as round a cylinder.
    static func grid(
        columns: Int,
        rows: Int,
        wrap: Bool = false,
        _ point: (Int, Int) -> Space3
    ) throws -> Mesh<Space3> {
        let across = wrap ? columns : columns + 1
        var vertices: [Space3] = []
        for j in 0...rows {
            for i in 0..<across { vertices.append(point(i, j)) }
        }
        func at(_ i: Int, _ j: Int) -> Index<Space3> { Index(Ordinal(UInt(j * across + i % across))) }
        var faces: [Polygon<Index<Space3>>] = []
        for j in 0..<rows {
            for i in 0..<columns {
                faces.append(try Polygon(vertices: [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)]))
            }
        }
        return try Mesh(vertices: vertices, faces: faces)
    }

    /// A rectangle `a` × `b` in the x–y plane with its lower left corner at the origin.
    static func rectangle(_ a: Double, _ b: Double, columns: Int, rows: Int) throws -> Mesh<Space3> {
        try grid(columns: columns, rows: rows) { i, j in
            Point(x: a * Double(i) / Double(columns), y: b * Double(j) / Double(rows), z: 0)
        }
    }

    /// The vertices on the outline of `rectangle(a, b, …)`.
    static func rim(_ mesh: Mesh<Space3>, _ a: Double, _ b: Double) -> [Int] {
        mesh.vertices.indices.filter { index in
            let p = mesh.vertices[index]
            return abs(p.x) < 1e-12 || abs(p.y) < 1e-12 || abs(p.x - a) < 1e-12 || abs(p.y - b) < 1e-12
        }
    }

    static func isotropic(density: Double, young: Double, poisson: Double) -> Acoustics.Plate.Material {
        .isotropic(
            density: Density.amount(density),
            young: Modulus.Young.pascals(young),
            poisson: Poisson.amount(poisson)
        )
    }

    static let steel = isotropic(density: 7_850, young: 210e9, poisson: 0.3)

    static func metres(_ value: Double) -> Measure<1, Space.Physical>.Value<Double> { Space.Physical.metres(value) }

    static func frequencies(_ modes: Acoustics.Shell.Modes) -> [Double] {
        modes.map { $0.frequency.underlying }
    }

    /// Each found frequency against the nearest expected one; the largest relative difference.
    static func worst(_ found: [Double], against expected: [Double]) -> Double {
        found.map { f in expected.map { abs(f / $0 - 1) }.min() ?? .infinity }.max() ?? .infinity
    }
}
