{ config, lib, pkgs, ... }:

let
  sci-python = pkgs.python3.withPackages (ps: with ps; [
      atlas
      numpy
      pandas
      olympus
      psycopg
      psycopg2
      pytorch
      seaborn
      scikit-learn
      sila2
      streamlit
  ]);
in
{
  environment.systemPackages = with pkgs; [
    sci-python
  ];
}
