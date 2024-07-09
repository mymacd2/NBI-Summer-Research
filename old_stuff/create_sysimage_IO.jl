using Pkg
using PackageCompiler

push!(LOAD_PATH, "/Users/millermacdonald/Desktop/Research_shit/NBI_Research")
using dsnbIOmodule

create_sysimage([:dsnbIOmodule], sysimage_path="sys_dsnbIOmodule.so")