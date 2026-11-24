public import Acoustics_Plate

extension Acoustics {
    /// A curved shell, such as a lute's bowl, a bell or a cymbal: the shape of its midsurface, how thick it is, what
    /// it is made of, the bars glued to it and how it is held. Its modes come from `Acoustics.Shell.Modes`.
    ///
    /// The midsurface is a mesh of quadrilaterals in space, coordinates in metres. Each face is a flat facet whose
    /// material has its first axis along `grain`, laid into the facet's plane: for a lute's bowl, the long axis of the
    /// body lies along every rib. A plate is a shell that happens to be flat.
    public struct Shell {
        /// The midsurface: quadrilateral faces over points in metres. A vertex no face uses is held.
        public var mesh: Mesh<Point<3, Double>>
        /// The thickness of every face not in `thicknesses`.
        public var thickness: Measure<1, Space.Physical>.Value<Double>
        /// The material of every face not in `materials`.
        public var material: Acoustics.Plate.Material
        /// The direction the materials' first axis follows, laid into each face's plane.
        public var grain: Vector<3, Double>
        /// Faces of another thickness, by their index in the mesh's faces.
        public var thicknesses: [Int: Measure<1, Space.Physical>.Value<Double>]
        /// Faces of another material, by their index in the mesh's faces.
        public var materials: [Int: Acoustics.Plate.Material]
        /// The bars glued to the surface.
        public var bars: [Bar]
        /// What is held at chosen vertices, by their index in the mesh's vertices.
        public var held: [Int: Hold]

        @inlinable
        public init(
            mesh: Mesh<Point<3, Double>>,
            thickness: Measure<1, Space.Physical>.Value<Double>,
            material: Acoustics.Plate.Material,
            grain: Vector<3, Double>,
            thicknesses: [Int: Measure<1, Space.Physical>.Value<Double>] = [:],
            materials: [Int: Acoustics.Plate.Material] = [:],
            bars: [Bar] = [],
            held: [Int: Hold] = [:]
        ) {
            self.mesh = mesh
            self.thickness = thickness
            self.material = material
            self.grain = grain
            self.thicknesses = thicknesses
            self.materials = materials
            self.bars = bars
            self.held = held
        }
    }
}

extension Acoustics.Shell: Hashable {}
extension Acoustics.Shell: Sendable {}
