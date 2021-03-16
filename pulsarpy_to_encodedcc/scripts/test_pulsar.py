import pulsarpy.models

model = getattr(pulsarpy.models, "Library")
rec = model.find_by({"id": 10721})
print(rec)
