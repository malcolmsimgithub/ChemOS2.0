{ config, lib, pkgs, ... }:
with pkgs;
{
  imports = [
    ./aiida.nix
  ];
}
