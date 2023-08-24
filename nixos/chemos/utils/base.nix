{ config, lib, pkgs, ... }:
with pkgs;
{
  imports = [
    ./nvim.nix
    ./overlay.nix
  ];
  environment.systemPackages = with pkgs; [
    graphviz
    openbabel
    killall
    snakemake
    tmux
    reptyr
    git
  ];
}
