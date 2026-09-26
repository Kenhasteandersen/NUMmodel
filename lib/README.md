This is where the compiled libraries are placed. They are not part of the
repository; build them from the fortran sources:

```sh
cmake -S . -B build
cmake --build build
cmake --install build
```

run from the root of the repository. This writes `libNUMmodel_matlab` (used by
the matlab code via `loadNUMmodelLibrary`) and `libNUMmodel_R` into this
directory, with the extension of the platform: `.dll` on windows, `.so` on
linux, `.dylib` on mac.

Compilation instructions are in the wiki.
