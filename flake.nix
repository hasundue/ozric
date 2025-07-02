{
  description = "An OZ equation solver using ceres-solver, cross-compiled with Zig";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";

    zig2nix = {
      url = "github:Cloudef/zig2nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    githooks-nix = {
      url = "github:cachix/git-hooks.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      nixpkgs,
      zig2nix,
      githooks-nix,
      treefmt-nix,
      ...
    }:
    let
      flake-utils = zig2nix.inputs.flake-utils;
    in
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        lib = nixpkgs.lib // {
          githooks-nix = githooks-nix.lib;
          treefmt-nix = treefmt-nix.lib;
        };

        # Zig flake helper
        # Check the flake.nix in zig2nix project for more options:
        # <https://github.com/Cloudef/zig2nix/blob/master/flake.nix>
        env = zig2nix.outputs.zig-env.${system} { };
      in
      with builtins;
      with env.pkgs.lib;
      rec {
        # Produces clean binaries meant to be ship'd outside of nix
        # nix build .#foreign
        packages.foreign = env.package {
          src = cleanSource ./.;

          # Packages required for compiling
          nativeBuildInputs = with env.pkgs; [
            ccache
            emscripten
          ];

          # Packages required for linking
          buildInputs = with env.pkgs; [ ];

          # Smaller binaries and avoids shipping glibc.
          zigPreferMusl = true;
        };

        # nix build .
        packages.default = packages.foreign.override (attrs: {
          # Prefer nix friendly settings.
          zigPreferMusl = false;

          # Executables required for runtime
          # These packages will be added to the PATH
          zigWrapperBins = with env.pkgs; [ ];

          # Libraries required for runtime
          # These packages will be added to the LD_LIBRARY_PATH
          zigWrapperLibs = attrs.buildInputs or [ ];
        });

        # For bundling with nix bundle for running outside of nix
        # example: https://github.com/ralismark/nix-appimage
        apps.bundle = {
          type = "app";
          program = "${packages.foreign}/bin/default";
        };

        # nix run .
        apps.default = env.app [ ] "zig build run -- \"$@\"";

        # nix run .#build
        apps.build = env.app [ ] "zig build \"$@\"";

        # nix run .#build-wasm
        apps.build-wasm =
          with pkgs;
          env.app [
            ccache
            emscripten
          ] "zig build --sysroot ${pkgs.emscripten}/upstream/emscripten wasm -- \"$@\"";

        # nix run .#test
        apps.test = env.app [ ] "zig build test -- \"$@\"";

        # nix run .#test-wasm
        apps.test-wasm = env.app [ pkgs.nodejs ] "node tests/test_wasm.js -- \"$@\"";

        # nix run .#docs
        apps.docs = env.app [ ] "zig build docs -- \"$@\"";

        # nix run .#zig2nix
        apps.zig2nix = env.app [ ] "zig2nix \"$@\"";

        # nix develop
        devShells.default =
          let
            treefmt = lib.treefmt-nix.mkWrapper pkgs {
              programs.nixfmt.enable = true;
              programs.zig.enable = true;
            };
            githooks = lib.githooks-nix.${system}.run {
              src = ./.;
              hooks.treefmt = {
                enable = true;
                package = treefmt;
              };
            };
          in
          env.mkShell {
            # Packages required for compiling, linking and running
            # Libraries added here will be automatically added to the LD_LIBRARY_PATH and PKG_CONFIG_PATH
            nativeBuildInputs =
              [ ]
              ++ packages.default.nativeBuildInputs
              ++ packages.default.buildInputs
              ++ packages.default.zigWrapperBins
              ++ packages.default.zigWrapperLibs;

            # Packages required for development
            packages = with pkgs; [
              clang-tools
              deno
              nodejs
              nil
              treefmt
              wasmtime
              zls
            ];
            shellHook = githooks.shellHook;
          };
      }
    );
}
