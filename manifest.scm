;; This "manifest" file can be passed to 'guix package -m' to reproduce
;; the content of your profile.  This is "symbolic": it only specifies
;; package names.  To reproduce the exact same profile, you also need to
;; capture the channels being used, as returned by "guix describe".
;; See the "Replicating Guix" section in the manual.

(specifications->manifest
  (list "r-seurat"
        "r-rstudioapi"
        "r-miniui"
        "r-styler"
        "r-rmarkdown"
        "r-stringr"
        "r-here"
        "r-matrix"
        "r-magrittr"
        "r-data-table"
        "r-dplyr"
        "r-deconvr"
        "libxml2"
        "r-rcurl"
        "r-guix-install"
        "lzlib"
        "r-tidyverse"
        "r-dyngen"
        "wget"
        "coreutils"
        "git"
        "bash"
        "findutils"
        "guix"
        "libxml2"
        "zlib"
        "curl"
        "openssl"
        "libxt"
        "cairo"
        "pkg-config"
        "cmake"
        "unzip"
        "tar"
        "pandoc"
        "which"
        "file"
        "diffutils"
        "gawk"
        "grep"
        "gzip"
        "gcc-toolchain@10"
        "gfortran-toolchain"
        "sed"
        "vim"
        "make"
        "python-pip"
        "python"
        "r"
        "r-guix-install"
        "glibc-locales"
        "nss-certs"))
