// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_xbio_extendr(void *dll);

void R_init_xbio(void *dll) {
    R_init_xbio_extendr(dll);
}
