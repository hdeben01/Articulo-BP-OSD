# 2026-08-28T17:11:55.264432200
import vitis

client = vitis.create_client()
client.set_workspace(path="hls")

comp = client.create_hls_component(name = "OSD_hls",cfg_file = ["hls_config.cfg"],template = "empty_hls_component")

comp = client.get_component(name="OSD_hls")
comp.run(operation="C_SIMULATION")

