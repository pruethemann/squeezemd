import openmm as mm
print([mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())])

platform = mm.Platform.getPlatformByName("CUDA")
print(platform)
print("CUDA version:", platform.getPropertyDefaultValue("CudaVersion"))
