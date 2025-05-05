# Bachelor Thesis: Single Energy Layer Proton Beam Range Modulator

This project strives to simulate a range modulator for single-energy-layer proton beams, which could be implemented in the research Software MatRAD and possibly in the clinical setting. The repository includes simulation scripts, modulator designs, CT/CST datasets, and tools for dose analysis.

---

## 🔧 Getting Started

Use the [`try_code.m`](try_code.m) script for examples on how to apply different modulators with pre-defined CTs and CSTs.

---

## 📁 Folder Structure

### `modulators/`
Contains all designed range modulators — both stable and experimental.

- `Box_modulator.m`  
  A stable, well-documented implementation with instructions on usage.

- `new_PTV_shaped_coned_inverse_good.m`  
  An experimental modulator that:
  - Takes the shape of the PTV.
  - Cuts a PTV-shaped region from the Box modulator.
  - Shapes cones to match the PTV outline.  
  ⚠️ Work in progress — parameter optimization needed for optimal results.

---

### `plane/`
Contains stable CT and CST functions, used for the already tested modulators.

---

### `data_analysis/`
Provides functions to:

- Analyze the dose produced by a given modulator.
- Compare dose distributions between different modulators.

These functions can be called at the end of a simulation script for data analysis.

---

### `plan/`
Includes functions that generate a single-energy-layer proton beam setup from the top of the CT.

---

## ✅ Notes

- Add all folders (`modulators`, `plane`, `data_analysis`, `plan`) to your MATLAB path.
- All components are modular — you can mix CTs, CSTs, and modulators as needed.
