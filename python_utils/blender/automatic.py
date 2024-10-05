import bpy
import sys



"""
mat = bpy.data.materials.new(name="potential")
mat.use_nodes = True
node_tree = mat.node_tree
nodes = node_tree.nodes
bsdf = nodes.get("Principled BSDF")
assert(bsdf)
vcol = nodes.new(type="ShaderNodeAttribute")
vcol.attribute_name = "Col"
node_tree.links.new(vcol.outputs[0], bsdf.inputs[0])
"""

#x = bpy.data.collections['meshes'].all_objects

# https://docs.blender.org/manual/en/latest/advanced/command_line/arguments.html

argv = sys.argv
argv = argv[argv.index("--") + 1:]  # get all args after "--"
print(argv)

#bpy.ops.wm.ply_import(filepath="/home/thomas/Downloads/dodecahedron.ply")
bpy.ops.import_mesh.ply(filepath=argv[0])
bpy.context.object.data.materials.append(bpy.data.materials['bubble'])

bpy.context.object.location[1] = 0.06
bpy.context.object.rotation_euler[0] = 1.53589

bpy.ops.object.modifier_add(type='BOOLEAN')
bpy.context.object.modifiers['Boolean'].object = bpy.data.objects['Cube']
bpy.ops.object.modifier_apply(modifier="Boolean")


#bpy.ops.transform.translate(value=(0,float(argv[0]),0))
