base_path <- "/space/"

rdn <- paste0(base_path, "Tyrol/ENMAP01-____L1C-DT0000001049_20220612T105735Z_028_V010303_20230922T131826Z_rdn_sub_small.img")
loc <- paste0(base_path, "Tyrol/ENMAP01-____L1C-DT0000001049_20220612T105735Z_028_V010303_20230922T131826Z_loc_sub_small.img")
obs <- paste0(base_path, "Tyrol/ENMAP01-____L1C-DT0000001049_20220612T105735Z_028_V010303_20230922T131826Z_obs_sub_small.img")
out_dir <- paste0(base_path, "test_lorema")
surface_path <- "/config_folder/surface_20240103_enmap_OG.json"

cmd <- sprintf(
  "docker exec isofit bash -c 'source ~/.bashrc && isofit apply_oe %s %s %s %s enmap --surface_path %s --emulator_base $(isofit path srtmnet --key file) --n_cores 10 --presolve'",
  rdn, loc, obs, out_dir, surface_path
)

# cmd <- sprintf(
#   "isofit apply_oe %s %s %s %s enmap --surface_path %s --emulator_base $(isofit path srtmnet --key file) --n_cores 10 --presolve",
#   rdn, loc, obs, out_dir, surface_path
# )
# 
# cmd <- sprintf(
#   "docker exec isofit bash -c 'source ~/.bashrc && isofit apply_oe %s %s %s %s enmap --surface_path %s --emulator_base /opt/isofit/srtmnet/sRTMnet_v120.h5 --n_cores 10 --presolve'",
#   rdn, loc, obs, out_dir, surface_path
# )
# 
# cmd <- sprintf(
#   "bash -c 'isofit apply_oe %s %s %s %s enmap --surface_path %s --emulator_base /opt/isofit/srtmnet/sRTMnet_v120.h5 --n_cores 10 --presolve'",
#   rdn, loc, obs, out_dir, surface_path
# )

# cmd <- sprintf(
#   " bash -c 'isofit apply_oe %s %s %s %s enmap --surface_path %s --emulator_base $(isofit path srtmnet --key file) --n_cores 10 --presolve'",
#   rdn, loc, obs, out_dir, surface_path
# )

print("start")
system(cmd)
print("emanuele")
