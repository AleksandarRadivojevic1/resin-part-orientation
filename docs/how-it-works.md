# How it works

This page explains what each command in `resin_orient.py` computes and why. The function names in brackets are where to look in the code.

## The setting

A resin printer of this kind builds the part hanging from a plate, one layer at a time, against a film at the bottom of the resin vat. After each layer the part is peeled off the film. Three things follow from that:

- A layer can only print if it touches something already printed. A region that starts with nothing above it on the plate side is an **island**. It cures against the film, stays there or drifts, and leaves a ragged edge.
- Peeling takes force in proportion to the **area of the layer**. Large layers late in a print pull on a part that is already long and heavy.
- A pocket whose opening faces the film is closed on all sides while its walls print: walls around it, a roof above it, the film below it. That is a **suction cup**. Each peel pulls a vacuum in it.

The tool describes a pose by one vector, `u`: the direction, in the part's own coordinates, that points away from the plate. Everything is then computed "along `u`". The previews draw the plate at the bottom, as a slicer does.

## Loading and mesh checks (`load_stl`, `load_3mf`, `cmd_info`)

Binary and ASCII STL are read directly. Identical vertices are merged, so the tool knows which triangles are neighbours. From a SatelLite `.goo.3mf` project it reads the mesh of one object; that mesh is already scaled and placed on the plate by the slicer.

`info` reports:

- size and volume;
- **open edges** (an edge with one triangle: a hole) and **non-manifold edges** (more than two triangles);
- the number of **shells**, meaning separate connected surfaces. Sculpted kits often consist of pieces that overlap without being merged. A shell with negative volume is an internal cavity or a piece that was exported inside out.

The two are told apart by where the shell lies (`fix_shells`): 24 points on its surface are tested against the other shells. If they are inside the rest of the part, the shell is a cavity and is kept. If not, it is a flipped piece and its triangles are re-wound. This matters because the slicer adds up surface crossings: left alone, a flipped piece counts as minus one, so it disappears and takes a bite out of whatever it overlaps. The fix is applied on every load, so `export` writes the corrected piece. A flipped piece buried entirely inside another one looks the same as a cavity and is left alone.

If the whole mesh is inside out, the tool flips it before analysing (`outward`).

## Sizes (`cmd_size`)

`size` takes every part of a kit, finds the box around all of them together and reads the height along the up axis (Z unless `--up` says otherwise). The scale for a size is the wanted height divided by that height: 150, 180 or 200 mm for S, M or L.

This only measures the figure if the files keep their assembled positions. Two checks guard against the common exceptions:

- if every part starts at the same height, the files are laid flat for printing, and the command stops and asks for `--native MM`;
- if all parts together are no taller than the tallest single part, it warns.

## The sweep (`cmd_sweep`, `island_heights`, `_persist`, `scores`)

The sweep scores 300 directions, spread evenly over the sphere by a Fibonacci spiral, plus the as-exported pose as the last row.

### Islands on the surface

For one direction, each vertex has a height above the plate. Imagine the resin level rising past the part from the plate side. Each local low point of the surface appears as a new, separate piece; as the level rises, pieces grow along the surface and merge.

The tool runs exactly that, with a union-find structure over the mesh edges sorted by height:

- A piece is **anchored** if it starts at the lowest point of the whole part, or at a low point whose surface faces away from the plate. The second case is the floor of a cavity: there is solid material beneath it, so nothing is hanging.
- Any other piece is an island. When it merges into an anchored piece, the height it has grown by then is how long it **floated**.

The sweep counts islands that float for more than 2 layers and for more than 5 layers. Shorter ones are bridged by normal curing.

This walk follows the surface, so it only sees connections within one shell. A piece on a second, overlapping shell never merges and is not counted. That is the reason for the `check` step.

### Overhang and support load

- **Down-facing area**: the area of triangles whose normal is within 45 degrees of pointing straight at the plate. This is the area that needs support.
- **Support load**: that area times its height above the plate plus the support lift. Tall supports under a large area cost resin and fail more often.
- **Height** of the part in this pose.

### Score

Each of the five numbers is divided by its median over all directions, so they are comparable, then weighted:

| Column | Weight |
|---|---|
| islands floating > 2 layers | 1 |
| islands floating > 5 layers | 3 |
| down-facing area | 1 |
| support load | 0.5 |
| height | 0.25 |

Lower is better. `rank --weights` changes the weights. The score narrows 300 directions to a handful; it does not make the final choice.

## The check (`cmd_check`, `slice_check`, `_crossings`)

`check` slices the solid in a given pose, at 0.3 mm pixels and one layer height by default.

When the poses come from a sweep, they are taken in score order but at least 25 degrees apart (`spaced`, `--spread`). Neighbouring sweep directions are about 12 degrees apart and score almost the same, so without this the slices go to several copies of one pose. The labels keep the sweep rank.

### Slicing

For each pixel, a vertical ray is tested against every triangle above that pixel (`_crossings`). Each crossing is recorded with the layer at which it happens and a sign: entering the solid or leaving it. Adding the signs up layer by layer gives, for every pixel of every layer, the number of surfaces the ray is inside.

A pixel is solid when that count is above zero. This is what makes overlapping shells work: where two pieces overlap the count is 2, which is still solid. An internal cavity brings the count back to 0, which is empty.

### Islands in the volume

In each layer, connected solid regions are labelled. A region is linked to the regions it overlaps in the layer before. Regions in the first layer are anchored. A region with nothing under it is an island, and it floats until it merges with an anchored region; a body that never merges floats to the top.

Because this works on the solid, not the surface, it counts islands correctly on parts made of several shells.

### Layer area

The solid area of every layer is the peel-force curve in the chart. `check` reports the largest layer, where it falls as a fraction of the height, and the largest layer in the last 20% of the height.

### Suction cups

In each layer the empty regions are labelled too. An empty region is **open** if:

- it touches the edge of the slice, meaning it connects to the resin around the part; or
- it is in the first layer, which opens to the gap under the part; or
- it overlaps an open region in the layer before.

Anything else is **sealed**: at the moment that layer prints, no path through the layers printed so far leads to open resin. Only printed layers count, because the layers still to come do not exist yet; a pocket that will later open at the top is still a cup while its walls print.

A pocket facing the plate overlaps the open space under it, so it is open from the start. A pocket facing the film starts under a solid roof and is sealed until its rim is reached. A hollow part with no drain holes is sealed throughout.

`check` reports the largest sealed area in any layer, the sealed volume, the height range and a grade. Up to 1 mm² is pixel noise and counts as no cup. A cup of up to 100 mm² and 0.3 ml is `minor`: that is the size of a key socket, and the Lessons in `SKILL.md` record sockets of 0.1 to 0.3 ml as tiny cups where a vent hole is optional. Anything larger is `MAJOR`. The two limits are `CUP_MINOR` in the code; they are a rule of thumb, not a measurement.

A pocket that faces the film is sealed at any tilt short of horizontal: its floor always prints before the lowest point of its rim. Tilting only makes the sealed part smaller. So a part with a socket at each end has a cup in nearly every pose, and the choice is between cup sizes.

`--cupscan N` looks for poses the score ranking would not offer. It slices the N best-scoring directions at 0.6 mm pixels and four layers per slice, which is enough to grade the cup, and then runs the full check on cup-free poses first and minor cups second, each in score order.

## Choosing a pose

The order of the rules in `SKILL.md`:

1. No `MAJOR` cup. No cup at all where a good pose allows it; a `minor` cup when it is a key socket and the pose is otherwise the best.
2. The largest layers should not come at the end of the print.
3. Glue faces toward the plate, tilted 20 to 45 degrees, so supports land where they will be hidden. Tips of scales and shingles away from the plate.
4. Then fewest long-floating islands, least visible down-facing area, least support.
5. Look at the preview. A pose that puts supports on a face or on fingers is rejected even when its numbers are better.

## Previews (`cmd_show`, `render`)

`show` rotates the part into the pose and draws it from several sides. Faces that need support are orange. Island starts found by the surface walk are circled in red when they are visible from that side.

The picture is drawn by scattering sample points over every triangle in proportion to its area on screen and keeping the nearest point per pixel. It is a preview, not a slicer view.

## Packing (`cmd_pack`, `footprint`)

For each part and each of 18 rotations about the vertical, the tool draws the footprint on a 0.5 mm grid, fills it, and grows it by the margin. Parts are placed largest first, each at the first free position found by testing the grown footprint against the plate with an FFT convolution. A beam search keeps the 24 tightest partial layouts at each step and the layout that uses the least plate length wins.

Since every footprint is grown by the margin, two parts end up about 1.5 to 2 times the margin apart, and a part is at least the margin plus 1 mm from the plate edge. `pack` warns when the tallest part plus the support lift exceeds the build height, and exits with an error when the parts do not fit.

## Export (`cmd_export`)

The part is scaled, rotated so `u` points up, turned about the vertical by the layout's angle, centred in X and Y, set down at z = 0 and written as binary STL. The slicer imports it at 100% with the tilt already applied.

## What the numbers do not tell you

- The sweep's island counts are approximate on multi-shell parts.
- `check` cannot see features thinner than its pixel. Halving the pixel (`--px 0.15`) roughly doubles the slicing time.
- Holes in the mesh let the ray count go wrong along a streak.
- No force is simulated. Layer area stands in for peel force; resin temperature, lift speed and film tension matter as much and are outside the tool.
