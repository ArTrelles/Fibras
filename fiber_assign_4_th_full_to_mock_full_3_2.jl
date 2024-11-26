include("../ModCoord.jl")
include("../ModKDTREE.jl")
include("../Modradec2xy.jl")

using .ModCoord # El punto (.) indica que es un módulo local
using .ModKDTREE 
using .Modradec2xy

using PyCall
using CSV
using DataFrames

using Base.Threads
using LinearAlgebra
using Random
#using StatsBase
using Dierckx

#para el kdtree
using NearestNeighbors

io = pyimport("io")
pandas = pyimport("pandas")
os = pyimport("os")
np = pyimport("numpy")
astropy = pyimport("astropy")
spa = pyimport("scipy.spatial")

#Carga ENV de python
ENV["PYTHON"] = "/home/ia/atrelles/miniconda3/envs/myenv/bin/python3"

#Carga tiles
load_tiles = pyimport("desimodel.io").load_tiles
tiles = load_tiles()

#convertir tiles a un df de Julia
tiles_df = tiles.to_pandas()
tiles_df.to_csv("tiles.csv")
tiles_df = CSV.File("tiles.csv") |> DataFrame

# Filtrar los tiles correspondientes al programa deseado, por ejemplo 'BRIGHT'
tiles_df = filter(row -> row.PROGRAM == "BRIGHT", tiles_df)

#Carga posiciones de las fibras

Table = astropy.table.Table
fiberpos_path = fiberpos.fits"

fiberpos_data = Table.read(fiberpos_path)

#println(fiberpos_data.colnames)

# Carga posiciones fibras como un df de pandas
fiberpos_df = fiberpos_data.to_pandas()

#crea un diccionario con fiberpos_df
fiberpos_dict = Dict(col => fiberpos_df[col].values for col in fiberpos_df.columns)

# convierte el diccionario a un DF de Julia
fiberpos_julia_df = DataFrame(fiberpos_dict)

fiberpos_petal = fiberpos_julia_df[:, 6]
fiberpos_x = fiberpos_julia_df[:, 12]
fiberpos_y = fiberpos_julia_df[:, 13]


#################funciones del desi model importadas a julia#######################


global _platescale = nothing

end



function get_radius_mm(theta::Vector{Float64})



function radec2xy_new(telra, teldec, ra, dec)


########Lectura de archivo de targets


# Definir el tamaño de la celda
cell_size = 1.5


using FITSIO

# Ruta al archivo FITS
input_fits = "/trabajo/atrelles/DESI_fibers_try/output/file.fits"

f = FITS(input_fits)

println("Número de HDUs: ", length(f))
#FITS(input_fits) do f
# Seleccionar directamente el segundo HDU, que contiene la tabla binaria
table_hdu = f[2]  # Cambia el índice según tu caso

if table_hdu isa FITSIO.TableHDU
        println("El HDU seleccionado contiene una tabla binaria.")

        # Convertir la tabla binaria a DataFrame
        sorted_df = DataFrame(table_hdu)

        # Mostrar información básica
        println("Columnas disponibles: ", names(sorted_df))
        println("Primeras filas:")
        println(first(sorted_df, 5))

        # Agregar nuevas columnas calculadas
        ra_targets = rad2deg.(sorted_df.ao1)
        dec_targets = rad2deg.(sorted_df.ao2)
        priorities = sorted_df.priot

        # Calcular cell_i y cell_j
        cell_size=1.5
        cell_i = floor.(Int, ra_targets / cell_size)
        cell_j = floor.(Int, dec_targets / cell_size)

        # Crear global_ids
        global_ids = collect(1:length(ra_targets))

        # Agregar columnas al DataFrame
        sorted_df.global_id = global_ids
        sorted_df.cell_i = cell_i
        sorted_df.cell_j = cell_j

        # Mostrar el DataFrame actualizado
        println("Número total de filas en `sorted_df`: ", nrow(sorted_df))
        println("Primeras filas de `sorted_df` actualizadas:")
        println(first(sorted_df, 5))
else
        println("El HDU seleccionado no es una tabla binaria.")
end


####################

# Fila con el `global_id` 4986
global_id_to_find = 4986

# Filtrar las filas con el `global_id` deseado
filtered_row = filter(row -> row.global_id == global_id_to_find, sorted_df)

# Imprimir la fila completa
if nrow(filtered_row) > 0
    println("Fila con global_id = $global_id_to_find encontrada:")
    println(filtered_row)
else
    println("No se encontró ninguna fila con global_id = $global_id_to_find.")
end

###################
# Cerrar el archivo FITS


ra_targets = convert(Vector{Float64}, ra_targets)
dec_targets = convert(Vector{Float64}, dec_targets)
priorities  = convert(Vector{Float64}, priorities)
cell_i = convert(Vector{Int}, cell_i)
cell_j = convert(Vector{Int}, cell_j)


# Agregar las columnas calculadas al DataFrame ordenado


#crea un diccionario Julia para cada celda de acuerdo a los targets que caen en la celda

cell_targets = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64, Float64, Int}}}()
for (ra, dec, priority, i, j, global_id) in zip(ra_targets, dec_targets, priorities, cell_i, cell_j, global_ids)
           cell = (i, j)
           if !haskey(cell_targets, cell)
               cell_targets[cell] = []
           end
           push!(cell_targets[cell], (ra, dec, priority, global_id))
       end
#println("Número de celdas en el diccionario: ", length(cell_targets))

kdtree = build_kdtree(fiberpos_x, fiberpos_y)
kdtree_py = build_kdtree_py(fiberpos_x, fiberpos_y)


cell_size=1.5 #grados
fb_len = length(fiberpos_x)
R_patrol=6.35


reserved_fibers_all_tiles = Vector{Tuple{Int64, Vector{Int64}}}()


#########Desacople de candidatos en cada tile
start_time_desacople = time()

# Diccionario para almacenar el tile_id asignado a cada global_id
target_to_tile = Dict{Int, Int}()


# Ahora target_to_tile contiene cada global_id con su tile_id asignado.
# Puedes proceder a la fase de asignación de fibras solo con los targets únicos.

end_time_desacople = time()
total_time_desacople = end_time_desacople - start_time_desacople
println("Tiempo total de ejecución: $total_time_desacople segundos")



####Primera asignacion (Necesita paralelizacion)
start_time = time()
##################################################
#using Base.Threads

# Inicializar las estructuras globales (vacías por ahora)
fiber_priorities_dict_global = Dict{Int64, Dict{Int64, Float64}}()
assignments_dict_global = Dict{Int64, Dict{Int64, Int64}}()
final_assignments_global = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64}}()

# Arreglo para almacenar los resultados temporales por tile (uno por hilo)
all_fiber_priorities = Vector{Dict{Int64, Dict{Int64, Float64}}}(undef, nthreads())
all_assignments = Vector{Dict{Int64, Dict{Int64, Int64}}}(undef, nthreads())
all_final_assignments = Vector{Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64}}}(undef, nthreads())


my_lock = ReentrantLock() 

# Paralelización de los tiles usando @threads
@threads for tile in eachrow(tiles_df)

    thread_id = Threads.threadid()

    # Crear diccionarios locales para este hilo
    local_fiber_priorities_dict = Dict{Int64, Dict{Int64, Float64}}()
    local_assignments_dict = Dict{Int64, Dict{Int64, Int64}}()
    local_final_assignments = Dict{Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64}}()

    tile_id = tile.TILEID
    tile_ra = tile.RA
    tile_dec = tile.DEC

    tile_ra_rad = deg2rad(tile_ra)
    tile_dec_rad = deg2rad(tile_dec)
 
    # Inicializar las prioridades de las fibras con -1
    #reserved_fibers = reserve_fibers(fiberpos_petal)
    
    # Agregar las fibras reservadas al arreglo de todos los mosaicos (tiles)
    #push!(reserved_fibers_all_tiles, (tile_id, reserved_fibers))


    reserved_fibers = get_reserved_fibers_from_vector(tile_id, reserved_fibers_all_tiles)
    if reserved_fibers === nothing
        println("No se encontraron fibras reservadas para tile_id=$tile_id")
        continue
    end


    # Encontrar targets candidatos
    cell_ra_idx = floor(Int, tile_ra / cell_size)
    cell_dec_idx = floor(Int, tile_dec / cell_size)
    cell_id = (cell_ra_idx, cell_dec_idx)

    candidate_target_ids = find_candidate_targets(tile_ra_rad, tile_dec_rad, cell_id, cell_targets, 0.0261799)

    # Proceso de asignación para cada target
    for global_id in candidate_target_ids
        #Solo procesar el global_id si está asignado a este tile en target_to_tile
        if haskey(target_to_tile, global_id) && target_to_tile[global_id] != tile_id
           continue  # Saltar si el global_id está asignado a otro tile
        end


        target_ra = ra_targets[global_id]
        target_dec = dec_targets[global_id]
        target_priority = priorities[global_id]

        target_ra_rad = deg2rad(target_ra)
        target_dec_rad = deg2rad(target_dec)

        target_x = 0.0
        target_y = 0.0


        try
           # Convertir RA y DEC del objetivo a coordenadas X, Y en el plano focal
         target_x, target_y = radec2xy_new(target_ra, target_dec, tile_ra, tile_dec)
        catch e
             # Si ocurre un error, se omite este objetivo
             #println("Error en radec2xy: ", e)
          continue
        end
        #println("   $global_id | $target_x | $target_y | $target_ra | $target_dec    ")


        #target_x, target_y = radec2xy_new(target_ra, target_dec, tile_ra, tile_dec)
        indices = inrange(kdtree, [target_x, target_y], R_patrol)

        for fiber_idx in indices
            if fiber_idx in reserved_fibers
                continue
            end
            if !haskey(local_assignments_dict, tile_id)
                local_assignments_dict[tile_id] = Dict(i => -1 for i in 1:fb_len)
            end
            if !haskey(local_fiber_priorities_dict, tile_id)
                local_fiber_priorities_dict[tile_id] = Dict{Int64, Float64}(i => -1 for i in 1:fb_len)
            end

            current_assignments = local_assignments_dict[tile_id]
            if current_assignments[fiber_idx] == -1
               fiber_x = fiberpos_x[fiber_idx]
               fiber_y = fiberpos_y[fiber_idx]
               distance = sqrt((target_x - fiber_x)^2 + (target_y - fiber_y)^2)
               if distance <= R_patrol
                                   local_final_assignments[(fiber_idx, tile_id)] = (tile_id, fiber_idx, global_id)
                  break
               end
            elseif target_priority > local_fiber_priorities_dict[tile_id][fiber_idx]
                #assigned_targets_count[thread_id] += 1 
                current_assignments[fiber_idx] = global_id
                                local_final_assignments[(fiber_idx, tile_id)] = (tile_id, fiber_idx, global_id)
                break
            end
        end
    end
    # Fusionar los resultados locales en los diccionarios globales dentro del bloqueo
    lock(my_lock) do
        merge!(fiber_priorities_dict_global, local_fiber_priorities_dict)
        merge!(assignments_dict_global, local_assignments_dict)
        merge!(final_assignments_global, local_final_assignments)
    end

    # Guardar los resultados locales en el arreglo global
    #all_fiber_priorities[thread_id] = local_fiber_priorities_dict
    #all_assignments[thread_id] = local_assignments_dict
    #all_final_assignments[thread_id] = local_final_assignments
end

# Sumar los contadores de todos los hilos
###total_assigned_targets = sum(assigned_targets_count)

###println("Total de targets asignados: $total_assigned_targets")
#=
#println("Total de targets asignados: $assigned_targets_count")
# Al final de la paralelización, fusionar los diccionarios locales en los globales
for local_fiber_priorities_dict in all_fiber_priorities
    merge!(fiber_priorities_dict_global, local_fiber_priorities_dict)
end

for local_assignments_dict in all_assignments
    merge!(assignments_dict_global, local_assignments_dict)
end

for local_final_assignments in all_final_assignments
    merge!(final_assignments_global, local_final_assignments)
end
=#
all_targets= Set(global_ids)

# Extraer la lista final de asignaciones
final_assignments_list = collect(values(final_assignments_global))
unassigned_targets = setdiff(all_targets, Set(t[3] for t in values(final_assignments_list)))
println("Número de targets no asignados: ",length(unassigned_targets))

assigned_global_ids = Set(t[3] for t in values(final_assignments_global))
println("Número de targets únicos asignados: ", length(assigned_global_ids))
##################################################
end_time = time()
total_time = end_time - start_time
println("Tiempo total final de ejecución: $total_time segundos")



################## TRABAJAR TARGETS NO ASIGNADOS
already_assigned_targets = Set{Int64}()

# Diccionario para almacenar los targets no asignados por celda
cell_targets_unassigned = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64, Float64, Int}}}()

# Recorrer solo los targets que están en unassigned_targets





for tile in eachrow(tiles_df)

    tile_id = tile.TILEID
    tile_ra = tile.RA
    tile_dec = tile.DEC

    # Convertir RA y DEC del tile a radianes
    tile_ra_rad = deg2rad(tile_ra)
    tile_dec_rad = deg2rad(tile_dec)

    # Encontrar la celda del tile actual
    cell_ra_idx = floor(Int, tile_ra / cell_size)
    cell_dec_idx = floor(Int, tile_dec / cell_size)
    cell_id =(cell_ra_idx,cell_dec_idx)

    # Verificar si el tile ya tiene asignaciones, si no, inicializarlas
    if !haskey(assignments_dict_global, tile_id)
         assignments_dict_global[tile_id] = Dict(i => -1 for i in 1:fb_len)
    end
    if !haskey(fiber_priorities_dict_global, tile_id)
            fiber_priorities_dict_global[tile_id] = Dict{Int64, Float64}(i => -1 for i in 1:fb_len)
    end

    current_assignments = assignments_dict_global[tile_id]   
 

    ## Buscar las fibras reservadas para el tile_id en reserved_fibers_all_tiles
    reserved_fibers = get_reserved_fibers_from_vector(tile_id, reserved_fibers_all_tiles)
    if reserved_fibers === nothing
        println("No se encontraron fibras reservadas para tile_id=$tile_id")
        continue
    end

    # Encontrar los targets candidatos para este tile (por ejemplo, usando la celda del tile)
    candidate_target_ids = find_candidate_targets(tile_ra_rad, tile_dec_rad, cell_id, cell_targets_unassigned, 0.0261799)

    for global_id in candidate_target_ids
       if global_id in already_assigned_targets
            continue  # Saltar objetivos que ya han sido asignados
          end

    # Recorrer los targets no asignados para intentar asignarlos a este tile
        # Obtener las coordenadas RA y DEC del target
        target_ra = ra_targets[global_id]
        target_dec = dec_targets[global_id]
        target_priority = priorities[global_id]

        target_x, target_y = 0.0, 0.0

        # Calcular las coordenadas X, Y del target en el plano focal
        try
            target_x, target_y = radec2xy_new(target_ra, target_dec, tile_ra, tile_dec)
        catch e
            continue  # Saltar si hay error en la conversión de coordenadas
        end

        # Buscar fibras dentro del radio de patrullaje usando el k-d tree de Python
        #indices = find_indices_py(kdtree_py, target_x, target_y, R_patrol)
        indices = inrange(kdtree, [target_x, target_y], R_patrol)


        # Intentar asignar el target no asignado a una fibra disponible
        for fiber_idx in indices
            if fiber_idx in reserved_fibers
                continue  # Saltar si la fibra está reservada
            end

            # Si la fibra está libre, asignar el target
            if current_assignments[fiber_idx] == -1
                fiber_x = fiberpos_x[fiber_idx]
                fiber_y = fiberpos_y[fiber_idx]
                distance = sqrt((target_x - fiber_x)^2 + (target_y - fiber_y)^2)

                if distance <= R_patrol
                                       final_assignments_global[(fiber_idx, tile_id)] = (tile_id, fiber_idx, global_id)
                    break  # Romper el bucle una vez que el target haya sido asignado
                end
                elseif target_priority > fiber_priorities_dict_global[tile_id][fiber_idx]
                  #if(global_id == 43954)
                    # println("entra el target 43954 a elseif una vez")
                  #end
                  displaced_global_id = current_assignments[fiber_idx] 
                                    push!(already_assigned_targets, global_id)
                  #println("$global_id global")
                  break
                end
        end
    end
end

# Recopilar los resultados finales en una lista
final_assignments_list = collect(values(final_assignments_global))
unassigned_targets_after_second_pass = setdiff(all_targets, Set(t[3] for t in final_assignments_list))
println("Número de targets no asignados después de la segunda fase: ", length(unassigned_targets_after_second_pass))

assigned_global_ids = Set(t[3] for t in values(final_assignments_global))
println("Número de targets únicos asignados al final de todo: ", length(assigned_global_ids))

##################

########Opcional
# Convertir final_assignments_list a un DataFrame

