public import Eigen

extension Acoustics.Shell {
    /// Why a shell's modes could not be found.
    public enum Error: Swift.Error {
        /// The face at this index is not a quadrilateral with positive area, or a thickness or material names no
        /// face.
        case face(Int)
        /// The bar at this index has no length, no section, or no side off the surface.
        case bar(Int)
        /// A held vertex that is not in the mesh.
        case vertex(Int)
        /// More modes asked for than the shell has unknowns.
        case count
        /// The eigenvalue solver failed.
        case eigen(Eigen.Error)
    }
}

extension Acoustics.Shell.Error: Hashable {}
extension Acoustics.Shell.Error: Sendable {}
