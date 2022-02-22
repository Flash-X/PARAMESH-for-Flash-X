Grid += \
  amr_build_bittree.o \
  amr_verify_bittree.o \
  amr_morton_order_bittree.o \
  amr_calculate_tree_data.o \
  amr_exchange_work_bflags.o \
  gr_btSortMortonBittree.o \
  gr_getIntCoords.o \
  gr_getNeighIntCoords.o


# The following GNU make special prevents that apparent dependencies
# on the file iso_c_binding.mod, which does usually not actually exist
# in the object directory but refers to the ISO_C_BINDING module known
# to FORTRAN compilers internally, trigger unnecessary recompilation
# of files that refer to the ISO_C_BINDING module.
.SECONDARY: iso_c_binding.mod
