{ config, lib, pkgs, ... }:

let
  mach-nix = import (builtins.fetchGit {
    url = "http://github.com/DavHau/mach-nix";
    rev = "8d903072c7b5426d90bc42a008242c76590af916";
  }) {
    pypiDataRev = "4cd2e59b40efb5848af9d64f95aff26fb23b6261";
    pypiDataSha256 = "sha256:19adgx3lbw8i1gczfqb4p5mk1gidrhzl9n95ij859a55hkvrjj66";
  };
  aiida-pkgs = mach-nix.mkPython {
    python = "python39";
    requirements = ''
      aiida-core
      aiida-orca
      aiida-gaussian
      aiida-shell
      pandas
      numpy
      scipy
    '';
  };

  verdi = pkgs.stdenv.mkDerivation {
    name = "verdi";
    buildInputs = [
      aiida-pkgs
    ];
    propagatedBuildInputs = [
      aiida-pkgs
    ];
    dontUnpack = true;
    installPhase = "install -Dm755 ${./verdi} $out/bin/verdi";
  };

in {
  environment.systemPackages = with pkgs; [ verdi ];
  services = {
    postgresql = {
      enable = true;
    };
    rabbitmq = {
      enable = true;
      configItems = {
        # Make it compatible with new AiiDA versions
        # https://github.com/aiidateam/aiida-core/wiki/RabbitMQ-version-to-use
        "consumer_timeout" = "3600000000";
      };
    };
  };
}
