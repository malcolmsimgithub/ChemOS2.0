{ config, lib, pkgs, ... }:

{
  nixpkgs.overlays = [(
    final: prev:
    {
      pythonPackagesOverlays = (prev.pythonPackagesOverlays or [ ]) ++ [
        (python-final: python-prev: {
          olympus = python-final.callPackage ./olympus.nix {};
          atlas = python-final.callPackage ./atlas.nix {};
          linear_operator = python-final.callPackage ./linear_operator.nix {};
          gpytorch = python-final.callPackage ./gpytorch.nix {};
          pyro-ppl = python-final.callPackage ./pyro-ppl.nix {};
          botorch = python-final.callPackage ./botorch.nix {};
          matter-golem = python-final.callPackage ./matter-golem.nix {};
          matter-chimera = python-final.callPackage ./matter-chimera.nix {};
          pyDOE = python-final.callPackage ./pyDOE.nix {};
          sobol-seq = python-final.callPackage ./sobol-seq.nix {};
        })
      ];
      python3 =
        let
          self = prev.python3.override {
          inherit self;
          packageOverrides = prev.lib.composeManyExtensions final.pythonPackagesOverlays;
        }; in
          self;
        python3Packages = final.python3.pkgs;
      }
    )
  ];
}
