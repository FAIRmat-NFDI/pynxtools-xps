# Examples for the Vamas coordinate system

Here, you can find an example for the coordinate system that is typically used in XPS and other surface science techniques. This is also the coordinate system which is defined by default in [NXxps](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html#nxxps).

The file `config_vms_cs_fixed.json` contains a hard-coded measurement configuration, which was translated to `NXxps` using the following command:

```console
user@box:~$ dataconverter config_vms_cs_fixed.json --reader xps --nxdl NXmpes --ignore-undocumented --output vms-cs-fixed.nxs 
```

The file `vms-cs.glb` contains a 3D representation (in [gltf/glb](https://en.wikipedia.org/wiki/GlTF) format) of the NXtransformation matrices in the NeXus file. It was created using the [nexus3d tool](https://github.com/domna/nexus3d) with the following command:

```console
user@box:~$ nexus3d vms-cs-fixed.nxs -fo vms-cs.glb --blender --left-handed
```

## Contact person in FAIRmat for these examples
Lukas Pielsticker