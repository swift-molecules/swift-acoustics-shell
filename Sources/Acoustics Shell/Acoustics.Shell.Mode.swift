public import Acoustics_Plate

extension Acoustics.Shell {
    /// One mode of a shell: how fast it rings, how each vertex moves, and the mass that moves with it.
    public struct Mode {
        /// The natural frequency.
        public let frequency: Frequency.Value<Double>
        /// How far each vertex of the mesh moves, in its order, scaled so that the longest is one; its largest
        /// component is positive. A held or unused vertex has zero.
        public let shape: [Vector<3, Double>]
        /// The modal mass at that scale: the mode's kinetic energy is one half this mass times the square of the
        /// speed of its farthest-moving vertex.
        public let mass: Mass.Value<Double>

        @inlinable
        public init(frequency: Frequency.Value<Double>, shape: [Vector<3, Double>], mass: Mass.Value<Double>) {
            self.frequency = frequency
            self.shape = shape
            self.mass = mass
        }
    }
}

extension Acoustics.Shell.Mode: Hashable {}
extension Acoustics.Shell.Mode: Sendable {}
