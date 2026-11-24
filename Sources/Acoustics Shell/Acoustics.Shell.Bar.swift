public import Acoustics_Plate

extension Acoustics.Shell {
    /// A straight bar of rectangular section glued to the surface, such as a soundboard's bar or a bridge: its own
    /// bending, stretching, twisting and shear about its centre line, which stands off the midsurface by half the
    /// shell's thickness and half its own height.
    ///
    /// It acts along the segment from `start` to `end`, laid into the plane of each face it crosses, over the part of
    /// it inside the face. A bar that runs along an edge between two faces counts half in each.
    public struct Bar {
        /// One end, on the midsurface.
        public var start: Point<3, Double>
        /// The other end, on the midsurface.
        public var end: Point<3, Double>
        /// The width across the surface.
        public var width: Measure<1, Space.Physical>.Value<Double>
        /// The height off the surface.
        public var height: Measure<1, Space.Physical>.Value<Double>
        /// The material, its first axis along the bar and its second across it in the surface: E_x stretches and
        /// bends it, G_xz shears it in its upright plane, G_xy twists it.
        public var material: Acoustics.Plate.Material
        /// Which side of the surface it stands on: any direction on that side.
        public var side: Vector<3, Double>

        @inlinable
        public init(
            from start: Point<3, Double>,
            to end: Point<3, Double>,
            width: Measure<1, Space.Physical>.Value<Double>,
            height: Measure<1, Space.Physical>.Value<Double>,
            material: Acoustics.Plate.Material,
            side: Vector<3, Double>
        ) {
            self.start = start
            self.end = end
            self.width = width
            self.height = height
            self.material = material
            self.side = side
        }
    }
}

extension Acoustics.Shell.Bar: Hashable {}
extension Acoustics.Shell.Bar: Sendable {}
