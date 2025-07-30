import h5py
import glob
import os

output_merged = "simdata.h5"
merged_file = h5py.File(output_merged, "w")

for filepath in sorted(glob.glob("./outputdata/simdata*.h5")):
    with h5py.File(filepath, "r") as f:
        v0 = f.attrs["v0"]
        group_name = f"v0_{v0:.2f}"
        data = f["data"][:]
        merged_file.create_dataset(group_name, data=data)
        merged_file[group_name].attrs["v0"] = v0
        print(f"Merged {filepath} as {group_name}")

merged_file.close()
print(f"✅ Merged data written to {output_merged}")

