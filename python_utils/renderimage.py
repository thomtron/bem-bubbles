import sys
import bpy

filename = sys.argv[1]
print(filename)

bpy.ops.import_mesh.ply(filepath="/home/thomas/Documents/ma-thesis/code/euler-results-pinned/pinned-beta/f=30e3_r=75e-6_p=6e3_beta=0.1_rem=0.12_epsilon=1e-2_b-nonlin-0.01/mesh-3000.ply")

bpy.context.object.active_material = bpy.data.materials.get("Gold")
bpy.ops.render.render(write_still = True)
bpy.ops.image.save_render(filepath="render-test.png")