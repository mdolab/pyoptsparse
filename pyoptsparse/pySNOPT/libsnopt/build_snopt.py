from argparse import ArgumentParser
import subprocess


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "source",
        type=str,
        help="The path to either the directory containing source files, or the path to the precompiled SNOPT library (libsnopt7.so, libsnopt7.dylib, or snopt7.dll) to link against",
    )
    parser.add_argument("--keep-build", action="store_true", help="Keep temporary build directory after building")
    parser.add_argument("--version", type=str, default="0.0.0", help="version of the snopt package")
    args = parser.parse_args()

    pip_commands = ["pip", "install", ".", f"-Csetup-args=-Dsourcedir={args.source}"]
    if args.keep_build:
        pip_commands += ["-Cbuild-dir=build"]
    subprocess.run(pip_commands, check=True)


if __name__ == "__main__":
    main()
