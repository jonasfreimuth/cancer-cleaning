;;; Run the following command to enter a development environment:
;;;
;;;  $ guix shell -m manifest.scm
;;;
;;; When the shell variable USE_GUIX_INFERIOR is set (to any value), a
;;; specific set of Guix channels will be used to build the
;;; environment.

(use-modules
 (guix profiles)
 (guix channels)
 (guix inferior)
 (ice-9 match))

(define channels
   (list (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "9e4632081ff31bf0d1715edd66f514614c6dc4bb")
        (introduction
          (make-channel-introduction
            "9edb3f66fd807b096b48283debdcddccfea34bad"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA"))))))

(define (lookup name)
  (specification->package name))

(define (lookup-inferior name)
  (define inferior
    (inferior-for-channels channels))
  (match (lookup-inferior-packages inferior name)
    ((first . rest) first)
    (_ (error
        (format #false "Could not find package `~a'.~%" name)))))

(define %packages
  (list "r-glmgampoi"
        "r-guix-install"
        "libxml2"
        "r-deseq2"
        "r-scuttle"
        "r-seurat"
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
        "r-rcurl"
        "lzlib"
        "r-tidyverse"
        "r-dyngen"
        "wget"
        "coreutils"
        "git"
        "bash"
        "findutils"
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
        "glibc-locales"
        "nss-certs"))

(packages->manifest
 (let ((how (if (getenv "USE_GUIX_INFERIOR")
                lookup-inferior lookup)))
   (map how (append %packages))))