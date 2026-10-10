#记录输入的模型
from materials.silicon import SILICON
from materials.base import Material
_Material = {
    "si" : SILICON,
    "silicon" : SILICON,
}
def get_material(name:str) -> Material:
    key = name.strip().lower()
    try:
        return _Material[key]
    except KeyError:
        raise ValueError(f"Unkonwn material:{name}")from None