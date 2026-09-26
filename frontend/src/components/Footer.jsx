import React from 'react';
import './footer.css';

const Footer = () => {
  return (
    <footer className="footer-strip">
      <div className="footer-strip-container">
        <p className="footer-text">
          © 2026 AQUADEX:{" "}
          <a href="mailto:ayaz.ahmd15@gmail.com" className="footer-link-item">Mail</a>
          {" "}|{" "}
          <a href="https://linkedin.com/in/ayazahmd" target="_blank" rel="noopener noreferrer" className="footer-link-item">LinkedIn</a>
          {" "}|{" "}
          <a href="https://github.com/Ayazxdev" target="_blank" rel="noopener noreferrer" className="footer-link-item">GitHub</a>
        </p>
        <p className="footer-subtext">
          Advancing environmental genomics through innovative AI solutions
        </p>
      </div>
    </footer>
  );
};

export default Footer;