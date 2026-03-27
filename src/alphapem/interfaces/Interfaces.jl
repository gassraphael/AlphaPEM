# -*- coding: utf-8 -*-

"""
    AlphaPEM.Interfaces

This module provides user-interface entry points for AlphaPEM, including the GUI
application launcher.
"""
module Interfaces

"""
    create_application()

Load the GUI script and launch the AlphaPEM graphical application.
"""
function create_application()
    # GUI.jl expects `Main.AlphaPEM` to exist when loaded.
    Base.include(Main, joinpath(@__DIR__, "GUI.jl"))
    return Main.create_application()
end

export create_application

end  # module Interfaces

