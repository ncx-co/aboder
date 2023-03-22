img <- magick::image_read("the_dude.jpeg")

hexSticker::sticker(
  img, 
  package="aboder",
  p_y = 0.2,
  p_color = "darkgreen",
  s_x = 1,
  s_y = 1.05,
  s_width = 1.5,
  s_height = 1.5,
  h_fill = "white",
  h_color = "darkgreen",
  filename = "man/figures/logo.png"
)
