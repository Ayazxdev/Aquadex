import React from 'react';
import './footer.css';

const Footer = () => {
  const currentYear = new Date().getFullYear();

  return (
    <footer className="footer-strip">
      <div className="footer-strip-container">
        <p className="footer-text">
          © {currentYear} AQUADEX by <strong>Hereka</strong>. All rights reserved.
        </p>
        <p className="footer-subtext">
          Advancing environmental genomics through innovative AI solutions
        </p>
      </div>
    </footer>
  );
};

export default Footer;