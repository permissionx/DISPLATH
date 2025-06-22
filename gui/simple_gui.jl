# Simple GUI for DISPLAΘ Ion Irradiation Simulator
# Requires: Pkg.add(["Gtk", "JSON3"])

using Gtk
using JSON3

# Include main simulation code
include("../src/main.jl")

# GUI Structure
mutable struct DISPLATH_GUI
    window::GtkWindow
    # Material Parameters
    lattice_constant_entry::GtkEntry
    box_size_x_entry::GtkEntry
    box_size_y_entry::GtkEntry
    box_size_z_entry::GtkEntry
    
    # Simulation Parameters  
    pmax_entry::GtkEntry
    vacancy_recover_entry::GtkEntry
    temperature_entry::GtkEntry
    debye_temp_entry::GtkEntry
    stop_energy_entry::GtkEntry
    
    # Ion Parameters
    ion_energy_entry::GtkEntry
    num_cascades_entry::GtkEntry
    
    # Material Selection
    material_combo::GtkComboBoxText
    ion_combo::GtkComboBoxText
    
    # Control Buttons
    run_button::GtkButton
    save_params_button::GtkButton
    load_params_button::GtkButton
    
    # Output
    output_text::GtkTextView
    progress_bar::GtkProgressBar
    
    # Status
    status_label::GtkLabel
end

function create_gui()
    # Create main window
    window = GtkWindow("DISPLAΘ - Ion Irradiation Simulator", 800, 600)
    
    # Create main vertical box
    main_vbox = GtkBox(:v, 10)
    set_gtk_property!(main_vbox, :margin, 10)
    
    # Title
    title_label = GtkLabel("DISPLAΘ - Binary Collision Approximation Simulator")
    set_gtk_property!(title_label, :markup, "<span size='large' weight='bold'>DISPLAΘ - Binary Collision Approximation Simulator</span>")
    
    # Create notebook for different parameter categories
    notebook = GtkNotebook()
    
    # === Material Parameters Tab ===
    material_grid = GtkGrid()
    set_gtk_property!(material_grid, :row_spacing, 5)
    set_gtk_property!(material_grid, :column_spacing, 10)
    set_gtk_property!(material_grid, :margin, 10)
    
    # Material selection
    material_combo = GtkComboBoxText()
    push!(material_combo, "SiC")
    push!(material_combo, "hBN") 
    push!(material_combo, "Graphene")
    push!(material_combo, "Custom")
    set_gtk_property!(material_combo, :active, 0)
    
    # Lattice parameters
    lattice_constant_entry = GtkEntry()
    set_gtk_property!(lattice_constant_entry, :text, "4.36")
    set_gtk_property!(lattice_constant_entry, :tooltip_text, "Lattice constant in Angstroms")
    
    # Box sizes
    box_size_x_entry = GtkEntry()
    box_size_y_entry = GtkEntry()
    box_size_z_entry = GtkEntry()
    set_gtk_property!(box_size_x_entry, :text, "50")
    set_gtk_property!(box_size_y_entry, :text, "50") 
    set_gtk_property!(box_size_z_entry, :text, "1400")
    set_gtk_property!(box_size_x_entry, :tooltip_text, "Number of unit cells in X direction")
    set_gtk_property!(box_size_y_entry, :tooltip_text, "Number of unit cells in Y direction")
    set_gtk_property!(box_size_z_entry, :tooltip_text, "Number of unit cells in Z direction")
    
    # Add to grid
    material_grid[1,1] = GtkLabel("Material:")
    material_grid[2,1] = material_combo
    material_grid[1,2] = GtkLabel("Lattice Constant (Å):")
    material_grid[2,2] = lattice_constant_entry
    material_grid[1,3] = GtkLabel("Box Size X:")
    material_grid[2,3] = box_size_x_entry
    material_grid[1,4] = GtkLabel("Box Size Y:")
    material_grid[2,4] = box_size_y_entry
    material_grid[1,5] = GtkLabel("Box Size Z:")
    material_grid[2,5] = box_size_z_entry
    
    push!(notebook, material_grid, "Material")
    
    # === Simulation Parameters Tab ===
    sim_grid = GtkGrid()
    set_gtk_property!(sim_grid, :row_spacing, 5)
    set_gtk_property!(sim_grid, :column_spacing, 10)
    set_gtk_property!(sim_grid, :margin, 10)
    
    pmax_entry = GtkEntry()
    vacancy_recover_entry = GtkEntry()
    temperature_entry = GtkEntry()
    debye_temp_entry = GtkEntry()
    stop_energy_entry = GtkEntry()
    
    set_gtk_property!(pmax_entry, :text, "4.0")
    set_gtk_property!(vacancy_recover_entry, :text, "4.0")
    set_gtk_property!(temperature_entry, :text, "300.0")
    set_gtk_property!(debye_temp_entry, :text, "490.0")
    set_gtk_property!(stop_energy_entry, :text, "10.0")
    
    set_gtk_property!(pmax_entry, :tooltip_text, "Maximum impact parameter")
    set_gtk_property!(vacancy_recover_entry, :tooltip_text, "Vacancy recovery distance")
    set_gtk_property!(temperature_entry, :tooltip_text, "Temperature in Kelvin")
    set_gtk_property!(debye_temp_entry, :tooltip_text, "Debye temperature in Kelvin")
    set_gtk_property!(stop_energy_entry, :tooltip_text, "Minimum energy in eV")
    
    sim_grid[1,1] = GtkLabel("pMax:")
    sim_grid[2,1] = pmax_entry
    sim_grid[1,2] = GtkLabel("Vacancy Recover Distance:")
    sim_grid[2,2] = vacancy_recover_entry
    sim_grid[1,3] = GtkLabel("Temperature (K):")
    sim_grid[2,3] = temperature_entry
    sim_grid[1,4] = GtkLabel("Debye Temperature (K):")
    sim_grid[2,4] = debye_temp_entry
    sim_grid[1,5] = GtkLabel("Stop Energy (eV):")
    sim_grid[2,5] = stop_energy_entry
    
    push!(notebook, sim_grid, "Simulation")
    
    # === Ion Parameters Tab ===
    ion_grid = GtkGrid()
    set_gtk_property!(ion_grid, :row_spacing, 5)
    set_gtk_property!(ion_grid, :column_spacing, 10)
    set_gtk_property!(ion_grid, :margin, 10)
    
    # Ion selection
    ion_combo = GtkComboBoxText()
    push!(ion_combo, "N")
    push!(ion_combo, "Ne")
    push!(ion_combo, "Ar")
    push!(ion_combo, "Kr")
    push!(ion_combo, "Xe")
    set_gtk_property!(ion_combo, :active, 0)
    
    ion_energy_entry = GtkEntry()
    num_cascades_entry = GtkEntry()
    
    set_gtk_property!(ion_energy_entry, :text, "10000.0")
    set_gtk_property!(num_cascades_entry, :text, "100")
    
    set_gtk_property!(ion_energy_entry, :tooltip_text, "Ion energy in eV")
    set_gtk_property!(num_cascades_entry, :tooltip_text, "Number of cascades to simulate")
    
    ion_grid[1,1] = GtkLabel("Ion Type:")
    ion_grid[2,1] = ion_combo
    ion_grid[1,2] = GtkLabel("Ion Energy (eV):")
    ion_grid[2,2] = ion_energy_entry
    ion_grid[1,3] = GtkLabel("Number of Cascades:")
    ion_grid[2,3] = num_cascades_entry
    
    push!(notebook, ion_grid, "Ion")
    
    # === Control Panel ===
    control_hbox = GtkBox(:h, 10)
    set_gtk_property!(control_hbox, :margin, 10)
    
    run_button = GtkButton("Run Simulation")
    save_params_button = GtkButton("Save Parameters")
    load_params_button = GtkButton("Load Parameters")
    
    set_gtk_property!(run_button, :name, "run-button")
    
    push!(control_hbox, run_button)
    push!(control_hbox, save_params_button)
    push!(control_hbox, load_params_button)
    
    # === Progress and Output ===
    progress_bar = GtkProgressBar()
    
    output_text = GtkTextView()
    set_gtk_property!(output_text, :editable, false)
    output_scroll = GtkScrolledWindow()
    set_gtk_property!(output_scroll, :min_content_height, 150)
    push!(output_scroll, output_text)
    
    # Status bar
    status_label = GtkLabel("Ready")
    set_gtk_property!(status_label, :halign, Gtk.GtkAlign.START)
    
    # Pack everything
    push!(main_vbox, title_label)
    push!(main_vbox, notebook)
    push!(main_vbox, control_hbox)
    push!(main_vbox, progress_bar)
    push!(main_vbox, GtkLabel("Output:"))
    push!(main_vbox, output_scroll)
    push!(main_vbox, status_label)
    
    push!(window, main_vbox)
    
    # Create GUI struct
    gui = DISPLATH_GUI(
        window,
        lattice_constant_entry, box_size_x_entry, box_size_y_entry, box_size_z_entry,
        pmax_entry, vacancy_recover_entry, temperature_entry, debye_temp_entry, stop_energy_entry,
        ion_energy_entry, num_cascades_entry,
        material_combo, ion_combo,
        run_button, save_params_button, load_params_button,
        output_text, progress_bar,
        status_label
    )
    
    # Connect signals
    setup_callbacks(gui)
    
    return gui
end

function setup_callbacks(gui::DISPLATH_GUI)
    # Run button callback
    signal_connect(gui.run_button, "clicked") do widget
        run_simulation(gui)
    end
    
    # Save parameters callback
    signal_connect(gui.save_params_button, "clicked") do widget
        save_parameters(gui)
    end
    
    # Load parameters callback
    signal_connect(gui.load_params_button, "clicked") do widget
        load_parameters(gui)
    end
    
    # Material combo callback
    signal_connect(gui.material_combo, "changed") do widget
        update_material_defaults(gui)
    end
    
    # Window close callback
    signal_connect(gui.window, "destroy") do widget
        Gtk.gtk_quit()
    end
end

function update_material_defaults(gui::DISPLATH_GUI)
    material = Gtk.bytestring(GAccessor.active_text(gui.material_combo))
    
    if material == "SiC"
        set_gtk_property!(gui.lattice_constant_entry, :text, "4.36")
        set_gtk_property!(gui.debye_temp_entry, :text, "490.0")
    elseif material == "hBN"
        set_gtk_property!(gui.lattice_constant_entry, :text, "2.51")
        set_gtk_property!(gui.debye_temp_entry, :text, "1000.0")
    elseif material == "Graphene"
        set_gtk_property!(gui.lattice_constant_entry, :text, "2.46")
        set_gtk_property!(gui.debye_temp_entry, :text, "2000.0")
    end
end

function get_parameters(gui::DISPLATH_GUI)
    return Dict(
        "lattice_constant" => parse(Float64, get_gtk_property(gui.lattice_constant_entry, :text, String)),
        "box_size_x" => parse(Int, get_gtk_property(gui.box_size_x_entry, :text, String)),
        "box_size_y" => parse(Int, get_gtk_property(gui.box_size_y_entry, :text, String)),
        "box_size_z" => parse(Int, get_gtk_property(gui.box_size_z_entry, :text, String)),
        "pmax" => parse(Float64, get_gtk_property(gui.pmax_entry, :text, String)),
        "vacancy_recover_distance" => parse(Float64, get_gtk_property(gui.vacancy_recover_entry, :text, String)),
        "temperature" => parse(Float64, get_gtk_property(gui.temperature_entry, :text, String)),
        "debye_temperature" => parse(Float64, get_gtk_property(gui.debye_temp_entry, :text, String)),
        "stop_energy" => parse(Float64, get_gtk_property(gui.stop_energy_entry, :text, String)),
        "ion_energy" => parse(Float64, get_gtk_property(gui.ion_energy_entry, :text, String)),
        "num_cascades" => parse(Int, get_gtk_property(gui.num_cascades_entry, :text, String)),
        "material" => Gtk.bytestring(GAccessor.active_text(gui.material_combo)),
        "ion_type" => Gtk.bytestring(GAccessor.active_text(gui.ion_combo))
    )
end

function log_message(gui::DISPLATH_GUI, message::String)
    buffer = get_gtk_property(gui.output_text, :buffer, GtkTextBuffer)
    insert!(buffer, message * "\n")
    
    # Scroll to end
    mark = get_gtk_property(buffer, :insert, GtkTextMark)
    GAccessor.scroll_mark_onscreen(gui.output_text, mark)
end

function update_status(gui::DISPLATH_GUI, status::String)
    set_gtk_property!(gui.status_label, :text, status)
end

function run_simulation(gui::DISPLATH_GUI)
    try
        update_status(gui, "Starting simulation...")
        log_message(gui, "=== Starting DISPLAΘ Simulation ===")
        
        # Get parameters
        params = get_parameters(gui)
        log_message(gui, "Parameters loaded: $(params["material"]) with $(params["ion_type"]) ions")
        
        # Set up simulation based on material
        if params["material"] == "SiC"
            # SiC setup
            a = params["lattice_constant"]
            primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
            boxSizes = [params["box_size_x"], params["box_size_y"], params["box_size_z"]]
            inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
            latticeRanges = [0 params["box_size_x"]; 0 params["box_size_y"]; 2 params["box_size_z"]-200]
            basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
            basisTypes = [1, 2]
            
            typeDict = Dict(
                1 => Element("Si", 43.0, 22.0),
                2 => Element("C", 40.0, 20.0),
                3 => Element(params["ion_type"], 1.0, 1.0)
            )
        else
            error("Material $(params["material"]) not yet implemented in GUI")
        end
        
        # Create parameters
        θτRepository = dirname(@__FILE__) * "/../thetatau_repository/"
        parameters = Parameters(
            primaryVectors, latticeRanges, basisTypes, basis,
            θτRepository, params["pmax"], params["vacancy_recover_distance"], 
            typeDict, 12345;
            temperature=params["temperature"], 
            DebyeTemperature=params["debye_temperature"],
            stopEnergy=params["stop_energy"]
        )
        
        log_message(gui, "Creating simulator...")
        simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)
        Save!(simulator)
        
        # Run cascades
        total_cascades = params["num_cascades"]
        vacancy_counts = Int[]
        
        for i in 1:total_cascades
            # Update progress
            progress = (i-1) / total_cascades
            set_gtk_property!(gui.progress_bar, :fraction, progress)
            update_status(gui, "Running cascade $i/$total_cascades")
            
            # Random ion position
            ionPosition = RandomPointInCircle(20.0) + [params["box_size_x"]*a/2, params["box_size_y"]*a/2, params["box_size_z"]*a*0.9]
            ion = Atom(3, ionPosition, parameters)
            SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
            SetEnergy!(ion, params["ion_energy"])
            push!(simulator, ion)
            
            # Run cascade
            Cascade!(ion, simulator)
            
            # Count vacancies
            vacancies, interstitials = DefectStatics(simulator)
            push!(vacancy_counts, length(vacancies))
            
            if i % 10 == 0
                log_message(gui, "Completed $i cascades, average vacancies: $(round(sum(vacancy_counts)/length(vacancy_counts), digits=2))")
            end
            
            # Process GUI events
            while Gtk.gtk_events_pending()
                Gtk.gtk_main_iteration()
            end
        end
        
        # Final results
        set_gtk_property!(gui.progress_bar, :fraction, 1.0)
        avg_vacancies = sum(vacancy_counts) / length(vacancy_counts)
        
        log_message(gui, "=== Simulation Complete ===")
        log_message(gui, "Total cascades: $total_cascades")
        log_message(gui, "Average vacancies per cascade: $(round(avg_vacancies, digits=2))")
        log_message(gui, "Vacancy count range: $(minimum(vacancy_counts)) - $(maximum(vacancy_counts))")
        
        # Save results
        results_file = "gui_results_$(now()).csv"
        open(results_file, "w") do f
            write(f, "cascade,vacancies\n")
            for (i, v) in enumerate(vacancy_counts)
                write(f, "$i,$v\n")
            end
        end
        log_message(gui, "Results saved to: $results_file")
        
        update_status(gui, "Simulation completed successfully")
        
    catch e
        log_message(gui, "ERROR: $e")
        update_status(gui, "Simulation failed")
        set_gtk_property!(gui.progress_bar, :fraction, 0.0)
    end
end

function save_parameters(gui::DISPLATH_GUI)
    try
        params = get_parameters(gui)
        filename = "displath_params_$(now()).json"
        open(filename, "w") do f
            JSON3.write(f, params)
        end
        log_message(gui, "Parameters saved to: $filename")
        update_status(gui, "Parameters saved")
    catch e
        log_message(gui, "Error saving parameters: $e")
    end
end

function load_parameters(gui::DISPLATH_GUI)
    try
        # Simple file dialog simulation - in practice, you'd use GtkFileChooser
        filename = "displath_params.json"
        if isfile(filename)
            params = JSON3.read(read(filename, String), Dict)
            
            # Update GUI fields
            set_gtk_property!(gui.lattice_constant_entry, :text, string(params["lattice_constant"]))
            set_gtk_property!(gui.box_size_x_entry, :text, string(params["box_size_x"]))
            set_gtk_property!(gui.box_size_y_entry, :text, string(params["box_size_y"]))
            set_gtk_property!(gui.box_size_z_entry, :text, string(params["box_size_z"]))
            set_gtk_property!(gui.pmax_entry, :text, string(params["pmax"]))
            set_gtk_property!(gui.vacancy_recover_entry, :text, string(params["vacancy_recover_distance"]))
            set_gtk_property!(gui.temperature_entry, :text, string(params["temperature"]))
            set_gtk_property!(gui.debye_temp_entry, :text, string(params["debye_temperature"]))
            set_gtk_property!(gui.stop_energy_entry, :text, string(params["stop_energy"]))
            set_gtk_property!(gui.ion_energy_entry, :text, string(params["ion_energy"]))
            set_gtk_property!(gui.num_cascades_entry, :text, string(params["num_cascades"]))
            
            log_message(gui, "Parameters loaded from: $filename")
            update_status(gui, "Parameters loaded")
        else
            log_message(gui, "File not found: $filename")
        end
    catch e
        log_message(gui, "Error loading parameters: $e")
    end
end

# Main function to run the GUI
function main()
    println("Starting DISPLAΘ GUI...")
    
    # Create and show GUI
    gui = create_gui()
    showall(gui.window)
    
    # Start GTK main loop
    if !Gtk.gtk_main_level()
        Gtk.gtk_main()
    end
end

# Auto-run if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end