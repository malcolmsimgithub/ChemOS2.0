{ config, lib, pkgs, ... }:
{
  environment.systemPackages = with pkgs; [
    postgresql
  ];
  services.postgresql = {
    enable = true;
    package = pkgs.postgresql;
    enableTCPIP = true;
    authentication = pkgs.lib.mkOverride 10 ''
      local all all trust
      host all all 127.0.0.1/32 trust
      host all all ::1/128 trust
    '';
    ensureUsers = [
    ];
    ensureDatabases = [
      "malcolm"
    ];
    initialScript = pkgs.writeText "backend-initScript" ''
      CREATE ROLE chemos WITH LOGIN PASSWORD 'chemos' CREATEDB;
      CREATE DATABASE chemos;
    '';
  };
}
