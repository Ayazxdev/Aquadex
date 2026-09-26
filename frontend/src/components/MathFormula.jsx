import React from "react";
import katex from "katex";
import "katex/dist/katex.min.css";

/**
 * High-performance, publication-grade LaTeX math formula renderer using KaTeX.
 * Supports inline ($...$) and block ($$...$$) mathematical rendering.
 */
export const MathFormula = ({ math, block = false, style = {} }) => {
  if (!math) return null;
  
  try {
    const html = katex.renderToString(math, {
      displayMode: block,
      throwOnError: false,
    });
    return (
      <span
        className={`math-formula ${block ? "math-block" : "math-inline"}`}
        style={{ display: block ? "block" : "inline-block", ...style }}
        dangerouslySetInnerHTML={{ __html: html }}
      />
    );
  } catch (err) {
    console.error("KaTeX render error:", err);
    return <span style={{ fontFamily: "monospace", ...style }}>{math}</span>;
  }
};

export default MathFormula;
