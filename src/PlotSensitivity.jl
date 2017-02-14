# includes -
using PyPlot
using PyCall
@pyimport matplotlib.patches as patches

# define my colors -
#color_2 = (1/255)*[255,238,225]
#color_1 = (1/255)*[255,99,71]

color_1 = (1/255)*[0,191,255]
color_2 = (1/255)*[240,248,255]

# How much of a gap between boxes?
epsilon = 0.2;

function draw_colormap(colorbar_col)

  # add the colorbar -
  for col_index = 1:1

    Z = colorbar_col

    # scale Z -
    Z = abs(Z)
    min_value = minimum(Z)
    max_value = maximum(Z)
    data_scaled = (Z-min_value)./(max_value - min_value)

    # how many patches per col?
    number_of_patches = length(data_scaled)
    ax = gca()
    for row_index = 1:number_of_patches

      origin_point = [(col_index - 1)+(col_index - 1)*epsilon+round(1.30*number_of_cols)+1,(row_index - 1)+(row_index - 1)*epsilon+1];

      # what color?
      fraction = (data_scaled[row_index]^2)/(0.15+data_scaled[row_index]^2)
      color_value = fraction*color_1+(1-fraction)*color_2


      # draw the square -
      ax[:add_patch](

                 patches.Rectangle(origin_point,   # (x,y)
                     1.0,          # width
                     1.0,          # height
                     facecolor=color_value,
                     edgecolor="grey",
                     linewidth=0.5,
                 )
             )


    end
  end
end

function draw_heatmap(data_slice)

  # main drawing loop -
  for col_index = 1:number_of_cols

    Z = data_slice[:,col_index]

    # scale Z -
    Z = abs(Z)
    min_value = minimum(Z)
    max_value = maximum(Z)
    if (max_value == min_value)

      data_scaled = zeros(15)

    else
      data_scaled = (Z-min_value)./(max_value - min_value)
    end

    @show size(data_scaled)

    # how many patches per col?
    number_of_patches = length(data_scaled)
    epsilon = 0.2;
    ax = gca()
    for row_index = 1:number_of_patches

      origin_point = [(col_index - 1)+(col_index - 1)*epsilon + 1,(row_index - 1)+(row_index - 1)*epsilon+1];

      # what color?
      fraction = (data_scaled[row_index]^2)/(0.15+data_scaled[row_index]^2)
      color_value = fraction*color_1+(1-fraction)*color_2

      # draw the square -
      ax[:add_patch](

                 patches.Rectangle(origin_point,   # (x,y)
                     1.0,          # width
                     1.0,          # height
                     facecolor=color_value,
                     edgecolor="grey",
                     linewidth=0.5,
                 )
             )


    end
  end
end

# Load the data -
data_array = readdlm("./sensitivity_results/avg_scaled_sensitivity_repressed.dat")

# what is the size of my array?
(number_of_rows,number_of_cols) = size(data_array)

# for col_index = 1:number_of_cols
#
#   # do we have all zeros in this col?
#   data_col = data_array[:,col_index]
#
#   idx_nz = find(data_col.== 0.0)
#   if (length(idx_nz) == number_of_rows)
#
#     for row_index = 1:number_of_rows
#       data_array[row_index,col_index] = 0.0001*rand()
#     end
#   end
# end

# add an extra col for colorbar -
colorbar_col = vec(transpose(linspace(0,1,number_of_rows)))

# draw the main plot -
draw_heatmap(data_array)

# draw the colorbar -
draw_colormap(colorbar_col)

# Make the axis look pretty -
axis("square")
axis("off")
savefig("./figs/Dynamic-lacI-Raw.pdf")
