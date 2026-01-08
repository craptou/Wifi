# IEEE 802.11a PHY Layer Implementation - AI Agent Instructions

## Project Overview
Implement the physical layer (PHY) of IEEE 802.11a WiFi standard, progressively adding complexity:
1. **Phase 1** (✅ Partial): OFDM over AWGN with ZF equalization
2. **Phase 2** (✅ Partial): OFDM over TGn-B multipath channel
3. **Phase 3** (In progress): Convolutional coding (7,1/2) + Viterbi decoding (hard/soft)
4. **Phase 4** (TODO): Dual-block interleaving per IEEE 802.11a spec
5. **Phase 5** (TODO): LDPC coding (rates 1/2, 2/3, 3/4, 5/6)
6. **Phase 6** (TODO): Channel estimation via 4 pilot subcarriers

**Architecture**: Parameter-centric functions with MATLAB Communications Toolbox. Focus on direct simulation without abstractions—every line must reflect IEEE 802.11a structure.

## Key Components & Data Flow

### 1. Parameter Management (`params_ofdm.m`)
**Single source of truth** per IEEE 802.11a Annex G:
- **OFDM**: FFT size $N=64$, useful subcarriers $N_u=48$, CP length $N_{cp}=16$
- **Modulation**: QPSK ($M=4, k=2$ bits/symbol). Phases 3+ allow 16-QAM, 64-QAM via `P.M`
- **Frequency**: $F_s = 20$ MHz (sampling), subcarrier spacing $\Delta f = 312.5$ kHz
- **TGn-B channel**: 9 multipath delays `tau_ns=[0 10 20 ...]` (ns), powers `p_dB=[0 -5.4 -10.8 ...]` (dB)
- **SNR sweep**: `EbN0_dB = 0:2:30`; each point runs independent SNR simulation
- **Scrambler seed**: `seed_canal=3` for reproducibility (used by `tgnb_channel.m` and future pilot generation)

### 2. Channel Model (`tgnb_channel.m`)
Generates **stationary Rayleigh-fading TGn-B realization**:
- Input: `P` (params), `seed` for RNG seeding
- Convert delays: $\tau_{samples}[k] = \text{round}(\tau_{ns}[k] / T_s)$ where $T_s = 1/F_s$
- Complex Gaussian coefficients: $\beta_k \sim \mathcal{CN}(0, p_{lin}[k])$ for each tap
- Output: temporal impulse response `h` (column vector), frequency response `H = FFT(h, Nfft)`, `Heff` (used post-equalization)
- **Pattern**: Each time function is called with same seed, produces identical `h`—critical for Phase 6 (pilot-based estimation)

### 3. OFDM Tx/Rx Simulation (`simulate_ofdm_tgnb.m`)
**Core processing pipeline** (must handle both coded/uncoded flows):

**Transmit side**:
1. Bit generation: `Nb = Nu * nb_symb_ofdm * k` (or doubled if coded)
2. **(Phase 3+)** Optional convolutional coding: `bits_coded = convenc(bits, trellis)` at rate 1/2
3. **(Phase 4+)** Optional interleaving: apply dual-block permutation (see below)
4. Modulation: `qammod(bits, M, 'InputType','bit', 'UnitAveragePower',true)` → constellation $X[m]$
5. Reshape to OFDM frame: $(N_{fft} \times N_{symb})$ matrix; place symbols on $N_u$ data subcarriers, zeros elsewhere
6. **(Phase 6+)** Pilot insertion: 4 fixed subcarriers with sequences from scrambler
7. IFFT + CP: `xt = ifft(Xf); xt_cp = [xt(end-Ncp+1:end,:); xt]` → time-domain vector

**Channel**:
- Convolution: `y_canal = conv(x, h); y_canal = y_canal(1:length(x))`
- AWGN: Noise power $P_{ne} = \frac{N}{N_u} \times \frac{P_{tx}}{2 \cdot k} \times \frac{1}{E_b/N_0}$ where $P_{tx}$ is avg TX power
  - This normalizes $E_b/N_0$ correctly across modulation orders and coded/uncoded comparisons
  
**Receive side**:
1. Reshape: `Ybloc = reshape(y, Nfft+Ncp, nb_symb)`
2. Remove CP: `yt = Ybloc(Ncp+1:end,:)`
3. FFT: `Yf = fft(yt)`
4. ZF equalization (Phase 2+): `Xhat_f = Yf ./ Heff` per subcarrier
   - MMSE available as alternative: $W_{mmse}[k] = H^*[k] / (|H[k]|^2 + \text{snr}^{-1})$
5. Extract data subcarriers: `z = Xhat_f(idx_data,:); z = z(:)` → flatten to vector
6. **(Phase 6+)** Channel estimation via pilots (before equalization step)
7. Soft demodulation (Phase 5+): `qamdemod(z, M, 'OutputType','llr', 'NoiseVariance', noiseVar)` → LLRs
   - For OFDM: $\tilde{\sigma}_k^2 = \sigma^2 / |H[k]|^2$ (effective noise variance per subcarrier post-equalization)
8. **(Phase 3+)** Hard demod or hard decisions: `qamdemod(z, M, 'OutputType','bit')` → binary
9. **(Phase 3+)** Viterbi decoding: `vitdec(bits_received, trellis, tb_depth, 'trunc', 'hard'/'soft')`
   - Traceback depth: $5 \times (K-1)$ where $K=3$ for (3,1/2) code, $K=7$ for (7,1/2)
10. TEB calculation: `TEB(ii) = mean(bits_est ~= bits_info)`

### 4. Interleaving (Phase 4) — `interleave_ofdm.m` (TODO)
**Dual-block interleaving per IEEE 802.11a § 2.3**:

**First interleaver** (frequency spreading):
$$k_s = \frac{N_{cbps}}{16} (k_e \bmod 16) + \lfloor k_e / 16 \rfloor, \quad k_e = 0,\ldots,N_{cbps}-1$$
where $N_{cbps} = N_u \times k$ (coded bits per OFDM symbol)

**Second interleaver** (within-symbol permutation):
$$k_s = s \lfloor k_e / s \rfloor + (k_e + N_{cbps} - \lfloor 16 k_e / N_{cbps} \rfloor) \bmod s$$
where $s = \max(k/2, 1)$ (bits per subcarrier ÷ 2)

**Implementation**: Pre-compute permutation indices; apply at encoder output before modulation. De-interleave at receiver post-demod, pre-decoder.

### 5. LDPC Coding (Phase 5) — `ldpc_encode.m`, `ldpc_decode.m` (TODO)
**Quasi-cyclic LDPC codes** per IEEE 802.11n:
- **Configuration**: Parity-check matrix $H$ defined by expansion factor $Z$ and permutation matrix $P$ (block circulant)
- **Rates**: 1/2, 2/3, 3/4, 5/6; block lengths 648, 1296, 1944
- **Example (R=1/2, N=648, Z=27)**: Use provided permutation matrix; build via `H = ldpcQuasiCyclicMatrix(blockSize, P)`
- **Encoding**: `cfgLDPCEnc = ldpcEncoderConfig(H); codeword = ldpcEncode(infoBits, cfgLDPCEnc)`
- **Decoding**: 
  - Config: `cfgLDPCDec = ldpcDecoderConfig(H, 'offset-min-sum'/'bp'/'layered-bp')`
  - Decode: `[decoded, niter, syndrome] = ldpcDecode(Lch, cfgLDPCDec, itermax, 'soft')`
  - Syndrome = 0 → convergence to valid codeword; else decoding failure
- **Key difference from Conv code**: Requires soft-decision **LLR** input; outputs LLRs if `DecisionType='soft'`

### 6. Pilot-Based Channel Estimation (Phase 6) — `estimate_channel_pilots.m` (TODO)
- **4 pilot subcarriers** at fixed positions (per 802.11a): subcarrier indices $\{7, 21, 43, 57\}$ (or similar)
- **Pilot data**: Generated via same LFSR as scrambler (but different init: all-ones)
- **Pilot symbols**: $\{+1, +1, +1, -1\}$ or $\{-1, -1, -1, +1\}$ per LFSR output
- **Estimation**: Linear interpolation or LS per pilot: $\hat{H}[i] = Y[i] / P[i]$ (received / known pilot)
- **Comparison**: Run simulation twice—once with perfect CSI, once with pilot-based estimate; compare TEB curves

### 7. Comparison Driver (`main_compare.m`)
- Generate **single** TGn-B realization via `seed_canal=3`
- Run `simulate_ofdm_tgnb(P, nb_symb, h, Heff, false)` → uncoded TEB
- Run same with `true` → coded TEB
- Plot curves: semilogy of TEB vs $E_b/N_0$ (dB)
- Annotate: Theoretical QPSK TEB (Phase 1), impact of multipath (Phase 2), coding gain (Phase 3+)

## Development Patterns & Workflow

### Phase Progression Strategy
1. **Phase 1**: Validate QPSK TEB vs. theoretical curve ($Q$-function for AWGN)—confirms modulation & noise model
2. **Phase 2**: Verify ZF equalization restores near-AWGN performance on multipath (magnitude response analysis)
3. **Phase 3**: Measure Conv(3,1/2) then Conv(7,1/2) coding gain; compare hard vs. soft Viterbi
4. **Phase 4**: Expect ~1-2 dB improvement from interleaving on bursty channel errors (breaks error correlation)
5. **Phase 5**: LDPC should outperform Conv codes at same rate; compare BP vs. min-sum algorithms
6. **Phase 6**: Pilot-only CSI estimation introduces ~1-3 dB penalty vs. perfect CSI at low SNR (channel variability)

### Parameter Modification Workflow
1. Edit `params_ofdm.m` (changes propagate downstream automatically)
2. If changing $F_e$, $\tau_{ns}$, or seed → regenerate channel: `[h, H, Heff] = tgnb_channel(P, seed)`
3. Re-run `simulate_ofdm_tgnb()` — MATLAB interprets (no compilation)
4. Phase 4+: Verify interleaver indices with small examples (4 bits, 2 subcarriers) before full run
5. Phase 5: Test LDPC with known test vectors first; validate syndrome = 0 on error-free channel

### Typical Extensions
- **Different modulation** (Phase 1→3): Replace `M=4` with 16 or 64; update `k=log2(M)` in `params_ofdm.m`
- **New channel** (Phase 2): Replace TGn-B in `tgnb_channel.m`, keep output format `[h, H, Heff]`
- **Conv code rate** (Phase 3): Use `poly2trellis(K, [octal_g1 octal_g2])`; update trellis, traceback depth
  - IEEE 802.11a standard: $K=7$, polynomials octal `[133 171]` (or punctured for rates 2/3, 3/4)
- **Soft-decision Viterbi** (Phase 3): Replace `vitdec(..., 'hard')` with `'soft'` and provide LLRs instead of bits
- **LDPC rates** (Phase 5): Pre-defined configs in IEEE 802.11n; load permutation matrix, call `ldpcEncoderConfig(H)`
- **Channel estimation** (Phase 6): Linear LS on pilots → interpolate to all subcarriers (FFT-based or simple linear)

### Testing Patterns (No formal test framework)
Validate subsystems manually:
- **Modulation**: Confirm $\mathbb{E}[|x[m]|^2] = 1$ after `qammod(..., 'UnitAveragePower', true)`
- **OFDM frame**: Check `Xf` shape $(N_{fft}, N_{symb})$; unused subcarriers are zero
- **Channel convolution**: Verify `length(y_canal) = length(x)` after truncation
- **Equalization**: Inspect `abs(Xhat_f)` distribution (should cluster near constellation points at high SNR)
- **TEB convergence**: TEB → 0 at high SNR for uncoded (confirms constellation recovery)
- **Coding gain**: Coded TEB should be $\approx 3$ dB lower than uncoded at same $E_b/N_0$ (Conv rate 1/2)
- **Interleaving**: Run with/without; expect noticeable gain on bursty errors
- **Pilot estimation**: CSI error ∝ channel SNR; visible as TEB floor at low SNR

## Important Conventions & Formulas

- **All signal vectors are column vectors** (rows = samples/subcarriers, columns = OFDM symbols)
- **Complex baseband throughout** (no carrier frequency; all processing at equivalent lowpass)
- **Power normalization**: `UnitAveragePower=true` ensures fair $E_b/N_0$ comparison across modulations
- **Noise model (AWGN per subcarrier)**: 
  $$P_{ne} = \frac{N_{fft}}{N_u} \times \frac{P_{tx}}{2k} \times \frac{1}{E_b/N_0}$$
  Divide $N_{fft}/N_u$ to account for spectral spreading; factor of $2k$ from complex symbols with $k$ bits each.
- **French comments preserved** (paramètres, canal, égalisateur, etc.) — maintain consistency in new code
- **No variable persistence**: Each SNR loop iteration is independent; no state between runs
- **Convolutional code standard**: Octal polynomials `[5 7]` for (3,1/2), `[133 171]` for IEEE 802.11a (7,1/2)
- **LLR computation (soft demod)**:
  $$L(b_i) = \log \frac{\sum_{x \in X_0^i} \exp(-|y/H - x|^2/\sigma_k^2)}{\sum_{x \in X_1^i} \exp(-|y/H - x|^2/\sigma_k^2)}$$
  where $\sigma_k^2 = \sigma^2 / |H[k]|^2$ (post-ZF effective noise variance)
- **Viterbi traceback depth**: Minimum $5(K-1)$ where $K$ is constraint length; use at least 10 for (7,1/2)

## File Naming & Output
- `.m` files: MATLAB source
- `.txt` files: Mirror `.m` (likely backups or version tracking)
- `results_compare_tgnb.mat`: Cached results (regenerate if parameters change)

## Common Pitfalls & Debugging
1. **Forgetting to reshape** after equalization from `(Nfft, nb_symb)` to vector for demod
2. **Noise power miscalculation**: Must account for `Nfft/Nu` scaling and `1/EbN0` relationship
3. **CP removal off-by-one**: Indices are `Ncp+1:end` to skip first `Ncp` samples
4. **Trellis definition mismatch**: Conv code uses octal `[5 7]` (simple) vs `[133 171]` (802.11a); Viterbi depth must match encoding constraint length $K$
5. **Interleaver index confusion**: First interleaver spreads bits across subcarriers; second permutes within each subcarrier. Pre-compute indices, don't recompute each symbol.
6. **LDPC syndrome non-zero**: If `syndrome ≠ 0`, decoder did not converge. Increase iterations or check LLR input range.
7. **Pilot position hardcoding**: Avoid magic numbers; define pilot indices in `params_ofdm.m` (e.g., `P.pilot_idx = [7 21 43 57]`)
8. **Channel estimation mismatch**: Perfect CSI in Phases 1-5 means `Heff = H`; Phase 6 replaces with estimated $\hat{H}$ only on data subcarriers.

## Reporting & Validation
- **Per phase comparison**: Always plot uncoded baseline, then coded variant on same figure with legend
- **Theoretical bounds**: AWGN QPSK (Phase 1) should match $Q$-function; multipath with ZF (Phase 2) should show 1-3 dB loss vs. AWGN
- **Coding gain**: Conv (1/2) typically 3-4 dB at TEB=1e-5; LDPC (1/2) typically 4-5 dB at same point
- **Interleaving impact**: Visible gain on bursty channels; minimal on stationary AWGN
- **Report figures**: All must have `title()`, `xlabel()`, `ylabel()`, `legend()` with proper units (e.g., "$E_b/N_0$ (dB)")
- **Code review**: Explain each function parameter, trace signal dimensions, justify choices vs. IEEE 802.11a standard
