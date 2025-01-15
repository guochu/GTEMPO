using Documenter, GTEMPO

makedocs(sitename="GTEMPO",
    pages=["index.md",
        "Manuals" => ["manuals/gtensor.md",
        "manuals/gmps.md"],
        "Examples" => ["examples/example1.md",
            "examples/example2.md"]])