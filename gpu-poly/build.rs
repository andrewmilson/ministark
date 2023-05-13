fn main() {
    // create a cfg alias for apple_silicon
    #[cfg(all(target_arch = "aarch64", target_os = "macos"))]
    println!("cargo:rustc-cfg=apple_silicon")
}
