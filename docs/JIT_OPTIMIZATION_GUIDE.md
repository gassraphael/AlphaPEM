# Guide d'Optimisation JIT avec PyCharm et Julia

## 📋 Vue d'ensemble

Ce guide explique comment optimiser les performances de compilation JIT (Just-In-Time) de Julia dans PyCharm et votre projet AlphaPEM.

---

## 1️⃣ Configuration de PyCharm pour Julia

### 1.1 Installation du Plugin Julia

1. **Ouvrir les Préférences** : `PyCharm` → `Preferences` (ou `Settings` sur Linux)
2. **Installer le Plugin Julia** :
   - Aller à `Plugins`
   - Chercher "Julia"
   - Installer le plugin officiel JetBrains Julia

### 1.2 Configuration de l'Interpréteur Julia

1. **Ajouter l'Interpréteur Julia** :
   - `Preferences` → `Languages & Frameworks` → `Julia`
   - Cliquer sur le bouton `...` à côté de "Julia executable"
   - Sélectionner votre exécutable Julia (généralement `/usr/bin/julia` ou `~/julia/bin/julia`)

2. **Vérifier la version** : 
   ```bash
   julia --version
   ```

### 1.3 Configurer les Options de Compilation JIT

Éditer le fichier `~/.julia/config/startup.jl` pour ajouter :

```julia
# Activer les optimisations
ENV["JULIA_THREADS"] = "auto"  # Utiliser tous les threads disponibles
ENV["JULIA_NUM_PRECOMPILE_TASKS"] = 4  # Nombre de tâches de précompilation

# Options de compilation JIT
if isdefined(Base.Experimental, :set_parameter_ex)
    Base.Experimental.set_parameter_ex(:JIT_COMPILER, :yes)
    Base.Experimental.set_parameter_ex(:JIT_OPTIMIZATION_LEVEL, :O2)
end
```

---

## 2️⃣ Stratégies de Compilation JIT dans Julia

### 2.1 Précompilation des Modules

Votre projet AlphaPEM utilise déjà la précompilation standard Julia. Pour améliorer les performances :

**Dans `src/AlphaPEM.jl`, assurez-vous que** :
- Tous les sous-modules sont correctement inclus
- Les dépendances lourdes (DifferentialEquations, Plots) sont lazy-loaded si possible

### 2.2 Type-Stabilité (Crucial pour JIT)

L'optimisation JIT dépend fortement de la **type-stabilité**. Exemple :

```julia
# ❌ PAS type-stable
function compute_unstable(x)
    if x > 0
        return x / 2  # Retourne Float64
    else
        return 0      # Retourne Int
    end
end

# ✅ Type-stable
function compute_stable(x)
    if x > 0
        return x / 2
    else
        return zero(typeof(x / 2))
    end
end
```

### 2.3 Utiliser `@time`, `@benchmark`, et `@profile`

```julia
using BenchmarkTools
using Profile

# Mesurer le temps JIT
@time result = run_simulation(cfg)  # Première exécution (compilation)
@time result = run_simulation(cfg)  # Deuxième exécution (code compilé)

# Benchmark détaillé
@benchmark run_simulation(cfg)

# Profiling
Profile.clear()
@profile for i in 1:10
    run_simulation(cfg)
end
Profile.print()
```

---

## 3️⃣ Précompilation Flexible (PackageCompiler.jl)

### 3.1 Installation

```julia
using Pkg
Pkg.add("PackageCompiler")
```

### 3.2 Créer une Image Système Précompilée

```julia
using PackageCompiler

create_sysimage(
    ["AlphaPEM", "Plots"],
    precompile_execution_file=joinpath(@__DIR__, "precompile_alphapem.jl"),
    sysimage_path="AlphaPEM.so",
    packages_as_names=true
)
```

### 3.3 Fichier de Précompilation (`precompile_alphapem.jl`)

```julia
# Simuler les cas d'usage courants pour la précompilation
import AlphaPEM.Config: SimulationConfig, PolarizationParams
import AlphaPEM.Application: run_simulation

# Cas 1: Polarization
current_params = PolarizationParams(
    delta_t_ini = 120.0 * 60.0,
    delta_i = 0.05e4,
    v_load = 0.01e4,
    delta_t_break = 15.0 * 60.0,
    i_max = 3.0e4
)

cfg = SimulationConfig(
    type_fuel_cell = :ZSW_GenStack,
    type_current = current_params,
    voltage_zone = :before_voltage_drop,
    type_auxiliary = :no_auxiliary,
    type_purge = :no_purge,
    type_display = :no_display,
    type_plot = :fixed
)

# Exécuter une itération pour déclencher la compilation
try
    run_simulation(cfg)
catch e
    @warn "Precompilation execution failed: $e"
end
```

### 3.4 Utiliser l'Image Système Précompilée

```bash
julia -J AlphaPEM.so examples/run_polarization.jl
```

---

## 4️⃣ Optimisations dans PyCharm

### 4.1 Configuration de Run dans PyCharm

**Créer une Julia Run Configuration** :

1. `Run` → `Edit Configurations`
2. Cliquer sur `+` et sélectionner `Julia`
3. Configuration suggérée :

```
Script path: /path/to/examples/run_polarization.jl
Arguments: (laisser vide pour les exemples)
Working directory: /path/to/AlphaPEM
Environment variables:
  JULIA_NUM_THREADS=auto
  JULIA_OPTIMIZE=2
```

### 4.2 Debugging Efficient en PyCharm

Pour déboguer sans pénalité de performance :

```julia
# Ajouter dans votre script
if haskey(ENV, "PYCHARM_DEBUGGING") || haskey(ENV, "JULIA_DEBUG")
    using Infiltrator  # Pour le debugging interactif
end
```

### 4.3 Utiliser les Breakpoints Intelligents

- Placer des breakpoints sur les lignes de traçage critiques
- Utilisar les "Conditional Breakpoints" pour éviter l'arrêt à chaque itération

---

## 5️⃣ Bonnes Pratiques pour AlphaPEM

### 5.1 Mesurer l'Impact du JIT

Modifiez vos scripts d'exemple pour inclure :

```julia
# Exemple amélioré: run_polarization.jl
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Config: SimulationConfig, PolarizationParams
using AlphaPEM.Application: run_simulation
using BenchmarkTools
using Logging

# Configuration
current_params = PolarizationParams(
    delta_t_ini = 120.0 * 60.0,
    delta_i = 0.05e4,
    v_load = 0.01e4,
    delta_t_break = 15.0 * 60.0,
    i_max = 3.0e4
)

cfg = SimulationConfig(
    type_fuel_cell = :ZSW_GenStack,
    type_current = current_params,
    voltage_zone = :before_voltage_drop,
    type_auxiliary = :no_auxiliary,
    type_purge = :no_purge,
    type_display = :no_display,
    type_plot = :fixed
)

# Première exécution: JIT compile
@info "First run (JIT compilation in progress)..."
@time run_simulation(cfg)

# Deuxième exécution: Code précompilé
@info "Second run (should be faster)..."
@time run_simulation(cfg)

# Benchmark statistique
@info "Benchmarking..."
stats = @benchmark run_simulation($cfg) samples=3
@info "Mean time: $(mean(stats.times) / 1e9) seconds"
```

### 5.2 Utiliser `@inbounds`, `@simd`, `@turbo` pour les Boucles Critiques

Identifiez les boucles dans `src/alphapem/core/models/` et optimisez :

```julia
# Exemple de boucle optimisée
@inbounds @simd for i in 1:length(data)
    result[i] = expensive_computation(data[i])
end
```

### 5.3 Lazy-Loading des Modules Lourds

Si Plots ne sont nécessaires que conditionnellement :

```julia
# Dans Application.jl
if config.type_display != :no_display
    import AlphaPEM.Utils.Display: display_results
    display_results(results)
end
```

---

## 6️⃣ Variables d'Environnement Importantes

```bash
# Pour les performances maximales
export JULIA_NUM_THREADS=auto
export JULIA_OPTIMIZE=2
export JULIA_NUM_PRECOMPILE_TASKS=4

# Pour le debugging
export JULIA_DEBUG=AlphaPEM

# Pour activer les mises en cache agressives
export JULIA_CODE_CACHE=yes
```

---

## 7️⃣ Diagnostic et Troubleshooting

### 7.1 Vérifier la Compilation JIT

```julia
# Voir les statistiques de compilation
@info "Compilation info:"
Base.print_precompile_statements()
```

### 7.2 Profiler la Compilation

```julia
using Profile
using ProfileSVG

prof = Profile.Allocs.@profile AlphaPEM.run_simulation(cfg)
ProfileSVG.savesvg("flame_graph.svg", prof)
```

### 7.3 Taille des Fichiers Précompilés

```bash
# Vérifier la taille du cache Julia
du -sh ~/.julia/compiled/v1.12/AlphaPEM/
```

---

## 8️⃣ Résumé des Actions à Faire

✅ **Immédiatement** :
1. Installer le plugin Julia dans PyCharm
2. Configurer l'interpréteur Julia
3. Ajouter `ENV["JULIA_NUM_THREADS"] = "auto"` à `~/.julia/config/startup.jl`

✅ **Court terme** :
1. Créer une Julia Run Configuration dans PyCharm
2. Mesurer les temps de compilation avec `@time`
3. Identifier les fonctions critiques et vérifier leur type-stabilité

✅ **Long terme** :
1. Utiliser PackageCompiler.jl pour créer une image système
2. Optimiser les boucles critiques avec `@simd` et `@inbounds`
3. Configurer le profiling continu avec GitHub Actions

---

## 📚 Ressources Supplémentaires

- [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
- [PackageCompiler.jl Documentation](https://julialang.github.io/PackageCompiler.jl/)
- [PyCharm Julia Plugin](https://plugins.jetbrains.com/plugin/10473-julia)
- [Julia Threading Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/)

---

**Dernière mise à jour** : 2026-04-05  
**Version Julia** : 1.12.5  
**Projet** : AlphaPEM v1.4.0

