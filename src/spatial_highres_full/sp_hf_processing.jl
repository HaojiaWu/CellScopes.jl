function convert_image_data(sp; layer_slot = "8_um")
    sp = set_default_layer(sp; layer_slot = layer_slot)
    sp = update_coordinates_hd(sp)
    return sp
end