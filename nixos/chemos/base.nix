{ config, lib, pkgs, ... }:
with pkgs;
{
  imports = [
    ./workflow_manager/base.nix
    ./experimental_planner/base.nix
    ./communication_layer/base.nix
    ./python.nix
    ./database/base.nix
    ./dft/base.nix
    ./utils/base.nix
    ./webgui/base.nix
  ];
}
