using Plotly

function plot_firing_rate(x, ys)
    data = []
    for y in ys
        trace = [
            "x"    => x,
            "y"    => y,
            "type" => "line"
        ]
        push!(data, trace)

data = [trace]
layout = [
  "showlegend" => false,
  "annotations" => [
    [
      "x" => 2,
      "y" => 5,
      "xref" => "x",
      "yref" => "y",
      "text" => "Annotation Text",
      "showarrow" => true,
      "arrowhead" => 7,
      "ax" => 0,
      "ay" => -40
    ]
  ]
]
response = Plotly.plot(data, ["layout" => layout, "filename" => "simple-annotation", "fileopt" => "overwrite"])
plot_url = response["url"]
