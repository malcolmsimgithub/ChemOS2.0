{ config, lib, pkgs, ... }:
with pkgs;
{
  imports = [
    ./postgresql.nix
  ];
}
