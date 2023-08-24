{ config, lib, pkgs, ... }:

{
  nixpkgs.overlays = [(
    final: prev:
    {
      pythonPackagesOverlays = (prev.pythonPackagesOverlays or [ ]) ++ [
        (python-final: python-prev: {
          singleton-decorator = python-final.callPackage ./singleton-decorator.nix {};
          pyunpack = python-final.callPackage ./pyunpack.nix {};
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
