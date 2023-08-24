{ config, lib, pkgs, ... }:
with pkgs;

{
  imports = [
    ./overlay.nix
  ];
}
