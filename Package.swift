// swift-tools-version: 6.4

import PackageDescription

let package = Package(
    name: "swift-acoustics-shell",
    platforms: [
        .macOS(.v27),
        .iOS(.v27),
        .tvOS(.v27),
        .watchOS(.v27),
        .visionOS(.v27),
    ],
    products: [
        .library(name: "Acoustics Shell", targets: ["Acoustics Shell"])
    ],
    dependencies: [
        .package(url: "https://github.com/swift-molecules/swift-acoustics.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-galerkin.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-sparse.git", branch: "main"),
        .package(url: "https://github.com/swift-molecules/swift-eigen.git", branch: "main"),
    ],
    targets: [
        .target(
            name: "Acoustics Shell",
            dependencies: [
                .product(name: "Acoustics", package: "swift-acoustics"),
                .product(name: "Galerkin", package: "swift-galerkin"),
                .product(name: "Sparse", package: "swift-sparse"),
                .product(name: "Eigen", package: "swift-eigen"),
            ]
        ),
        .testTarget(
            name: "Acoustics Shell Tests",
            dependencies: [
                .target(name: "Acoustics Shell")
            ]
        ),
    ],
    swiftLanguageModes: [.v6]
)

for target in package.targets where ![.system, .binary, .plugin, .macro].contains(target.type) {
    target.swiftSettings =
        (target.swiftSettings ?? []) + [
            .strictMemorySafety(),
            .enableUpcomingFeature("ExistentialAny"),
            .enableUpcomingFeature("InternalImportsByDefault"),
            .enableUpcomingFeature("MemberImportVisibility"),
            .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
            .enableUpcomingFeature("InferIsolatedConformances"),
            .enableExperimentalFeature("Lifetimes"),
            .treatAllWarnings(as: .error),
        ]
}
