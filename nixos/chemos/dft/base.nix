{ config, lib, pkgs, ... }:
with pkgs;
let
  multiwfn = pkgs.callPackage ./multiwfn.nix {};
in
{
  environment.systemPackages = with pkgs; [
    multiwfn
  ];
}
