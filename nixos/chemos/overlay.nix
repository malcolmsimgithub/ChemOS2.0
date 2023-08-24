{ config, lib, pkgs, ... }:

{
  nixpkgs.overlays = [(
    final: prev:
    {
      pythonPackagesOverlays = (prev.pythonPackagesOverlays or [ ]) ++ [
        (python-final: python-prev: {
          sila2 = python-final.callPackage ./communication_layer/sila2.nix {};
          olympus = python-final.callPackage ./experimental_planner/olympus.nix {};
          atlas = python-final.callPackage ./experimental_planner/atlas.nix {};
          linear_operator = python-final.callPackage ./experimental_planner/linear_operator.nix {};
          gpytorch = python-final.callPackage ./experimental_planner/gpytorch.nix {};
          pyro-ppl = python-final.callPackage ./experimental_planner/pyro-ppl.nix {};
          botorch = python-final.callPackage ./experimental_planner/botorch.nix {};
          matter-golem = python-final.callPackage ./experimental_planner/matter-golem.nix {};
          matter-chimera = python-final.callPackage ./experimental_planner/matter-chimera.nix {};
          pyDOE = python-final.callPackage ./experimental_planner/pyDOE.nix {};
          sobol-seq = python-final.callPackage ./experimental_planner/sobol-seq.nix {};
          streamlit = python-final.callPackage ./webgui/streamlit.nix {};
          singleton-decorator = python-final.callPackage ./utils/singleton-decorator.nix {};
          pyunpack = python-final.callPackage ./utils/pyunpack.nix {};
        })
      ];
      python3 =
        let
          self = prev.python310.override {
          inherit self;
          packageOverrides = prev.lib.composeManyExtensions final.pythonPackagesOverlays;
        }; in
          self;
        python3Packages = final.python3.pkgs;
      }
    )
  ];
}
