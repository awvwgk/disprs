wide_arithmetic = host_machine.cpu_family() in [
  'aarch64', 'loongarch64', 'mips64', 'ppc64', 'riscv64', 'wasm32', 'x86_64',
]
if not wide_arithmetic
  wide_arithmetic = meson.get_compiler('rust').compiles(
    '#[cfg(not(target_pointer_width = "64"))] compile_error!("not 64-bit"); fn main() {}',
    name: '64-bit pointer width')
endif
extra_args += ['--cfg', 'fast_arithmetic="' + (wide_arithmetic ? '64' : '32') + '"']