# dec406 hex decoder (bundled)

The observatory reads each frame's fields and the position a beacon reports with the hex decoder
of [dec406](https://github.com/moricef/dec406) by Fabrice Morel, MIT licence (see `LICENSE` here).

These files are copied unchanged from dec406 at commit
`b811499da3b239fd37f5744fa6674dab48a1b9a4` (2026-10-30): `src/dec406_hex.c`, `src/dec406.c`,
`src/dec406_v1g.c`, `src/dec406_v2g.c`, `src/display_utils.c` and the headers they include.
CMake builds them into `dec406_hex`, which the installer puts in `/opt/meosar-rx/bin`. The
receiver itself compiles none of this code.

To update: copy the same files from a newer dec406 and change the commit above.
