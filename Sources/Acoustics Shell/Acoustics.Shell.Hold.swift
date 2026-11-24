extension Acoustics.Shell {
    /// The motions held at a vertex, along and about the axes of space.
    public struct Hold: OptionSet {
        public let rawValue: UInt8

        @inlinable
        public init(rawValue: UInt8) {
            self.rawValue = rawValue
        }

        /// No motion along x.
        public static let x = Self(rawValue: 1 << 0)
        /// No motion along y.
        public static let y = Self(rawValue: 1 << 1)
        /// No motion along z.
        public static let z = Self(rawValue: 1 << 2)
        /// No turning about x.
        public static let turnX = Self(rawValue: 1 << 3)
        /// No turning about y.
        public static let turnY = Self(rawValue: 1 << 4)
        /// No turning about z.
        public static let turnZ = Self(rawValue: 1 << 5)
        /// The vertex does not move but may turn.
        public static let pinned: Self = [.x, .y, .z]
        /// The vertex neither moves nor turns.
        public static let clamped: Self = [.x, .y, .z, .turnX, .turnY, .turnZ]
    }
}

extension Acoustics.Shell.Hold: Hashable {}
extension Acoustics.Shell.Hold: Sendable {}

extension Acoustics.Shell.Hold {
    /// The unknowns held, as components 0–5: the motions along x, y and z, then the turns about them.
    var components: Set<Int> {
        Set((0..<6).filter { rawValue & (1 << $0) != 0 })
    }
}
