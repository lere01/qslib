fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    match qslib_cli::run(&args) {
        Ok(output) => println!("{output}"),
        Err(error) => {
            eprintln!("qslib: {error}");
            std::process::exit(2);
        }
    }
}
