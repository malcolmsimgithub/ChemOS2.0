{ config, lib, pkgs, ... }:

{
  programs = {
    neovim = {
      enable = true;
      defaultEditor = true;
      vimAlias = true;
      viAlias = true;
      configure = {
        customRC = ''
          set nu
        '';
        packages.myVimPackage = with pkgs.vimPlugins; {
        start = [
          fugitive
          vim-nix
          vim-airline
          vim-airline-themes
          ale
        ];
        };
      };
    };
  };
}
